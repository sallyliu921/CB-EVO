#pragma once
#include <vector>
#include <string>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <unordered_map>
#include <random>
#include <fmt/core.h>
#include "cmd.hpp"
#include "utils.hpp"

namespace fs = std::filesystem;

// Global settings and paths
fs::path tmp_dir_path{"./tmp"};

// Default optimization commands
std::string add_aig_opt_cmd = "strash; ifraig; dch -f; strash"; // Add once every 12 stages.
std::string add_fpga_opt_cmd = "strash; ifraig; scorr; dc2; strash; dch -f; if -K 6; mfs2; lutpack -S 1; strash"; // Add once every 12 stages.
std::string fpga_cmd_flowtune = "if -K 6; ";
std::string fpga_cmd = "if -a -K 6; ";
std::vector<std::string> vtr_cases = {"bfly", "dscg", "fir", "or1200","ode","syn2"};

// Global variables for command line arguments
std::string input_file;
int target = 0;
int seq_all = 0;
int exp_num = 0;
bool use_static_extraction = false; 
double time_limit;

// init node definitions
struct Oper{
    Cmd type;
    uint32_t arg0;
    bool f0 = false;
    Oper(Cmd t = Cmd::B, uint32_t a0 = 0): type(t), arg0(a0){}
};

// Function definitions
struct Node;
Cmd stringToCmd(const std::string& cmd);
Oper str2oper(const std::string& soper);
std::string oper2str(CmdManager& manager, const Oper &oper);
Oper rand_oper(std::vector<Cmd> types);
Node get_ori_node(const fs::path& i_path);
// Oper rand_oper_map(Cmd type);

// Helper function to set time limit based on AIG node count
inline void set_time_limit_by_aig_nodes(int aig_node_count) {
    if (aig_node_count < 1000) {
        time_limit = 60;
    } else if (aig_node_count < 10000) {
        time_limit = 300;
    } else {
        time_limit = 3600;
    }
}


/////////////////////// Node structure definitions ///////////////////////
struct Node{
    fs::path aig_addr;   
    int length = 0;
    std::vector<Oper> seq; // exclude fpga command

    long aig_node = -1, aig_lev = -1;
    long lut_area = -1, lut_lev = -1;

    // Calculate cost based on target type
    inline int calculate_cost() const {
        if (target == 0) {
            return aig_node;
        } else if (target == 1) {
            return aig_node * aig_lev;
        } else if (target == 2) {
            return lut_area;
        } else {  // target == 3
            return lut_area * lut_lev;
        }
    }

    // Static version for use without a Node instance
    static inline int calculate_cost(const Node& node) {
        return node.calculate_cost();
    }

    std::string get_seq_str(CmdManager& manager) const{
        std::string res;
        for(const auto &oper : seq){
            res += oper2str(manager, oper) + "; "; 
        }
        res += (target >= 2) ? ((exp_num == 1) ? fpga_cmd_flowtune : fpga_cmd) : "";
        return res;
    }

    std::unordered_map<Cmd, Node> get_short_term_reward(Node& node, CmdManager& manager){
        std::unordered_map<Cmd, Node> map_st_nodes;
        std::unordered_map<Cmd, std::string> map_cmd_st = constr_short_term_commands(node, manager);
        
        std::vector<std::pair<Cmd, std::string>> cmd_list(map_cmd_st.begin(), map_cmd_st.end());
        std::vector<std::pair<Cmd, Node>> results(cmd_list.size());
        fs::path st_save_path = node.aig_addr.parent_path();

        #pragma omp parallel for
        for (int i = 0; i < cmd_list.size(); ++i) {
            Node tmp_node;
            std::string output = executeAbcCommand(cmd_list[i].second);
            if (target == 0 || target == 1){
                aig_res ress = res_aig(remove_ansi_sequences(output));
                tmp_node.aig_node = ress.aig_node;
                tmp_node.aig_lev = ress.levs;
            } else {
                map_res ress = res_lut(remove_ansi_sequences(output));
                tmp_node.lut_area = ress.lut_nd;
                tmp_node.lut_lev = ress.lut_lev;
            }

            tmp_node.aig_addr = st_save_path / fmt::format("cb_st_{}.aiger", static_cast<int>(cmd_list[i].first));

            #pragma omp critical
            {
                results[i] = {cmd_list[i].first, tmp_node};
            }
        }

        for (const auto& result : results) {
            map_st_nodes[result.first] = result.second;
        }
        
        return map_st_nodes;
    }

    std::unordered_map<Cmd, std::vector<Node>> get_long_term_reward(Node& node, CmdManager& manager, int search_times, const std::vector<int>& search_depths)
    {
        int numArms = manager.getCurrentMapSize();
        auto cmdsMap = manager.getCurrentTypeToStr();

        std::unordered_map<Cmd, std::vector<std::string>> map_cmd_lt = constr_long_term_commands(node, manager, search_times, search_depths);
        
        fs::path lt_save_path = node.aig_addr.parent_path();
        int total_length=0;
        int groupSize = search_times * search_depths.size();

        std::unordered_map<Cmd, std::vector<Node>> map_lt_nodes;
        for (const auto& c : cmdsMap) {
            Cmd arm = c.first;

            std::vector<Node> nodes(search_depths.size()*search_times);  
            map_lt_nodes[arm] = nodes;
            
            #pragma omp parallel for
            for (int i = 0; i < search_depths.size()*search_times; i++){

                Node tmp_node;
                std::string output = executeAbcCommand(map_cmd_lt[arm][i]);
                if (target == 0 || target == 1){
                    aig_res ress = res_aig(remove_ansi_sequences(output));
                    tmp_node.aig_node = ress.aig_node;
                    tmp_node.aig_lev = ress.levs;
                } else {
                    map_res ress = res_lut(remove_ansi_sequences(output));
                    tmp_node.lut_area = ress.lut_nd;
                    tmp_node.lut_lev = ress.lut_lev;
                }

                #pragma omp critical
                {
                    map_lt_nodes[arm][i] = tmp_node;
                }
            }

            std::string rm_tmp_aig = fmt::format("rm  {}/cb_lt_{}*", lt_save_path.string(), arm);
            system(rm_tmp_aig.c_str());
        }   

        return map_lt_nodes;
    }
    
    std::unordered_map<Cmd, std::string> constr_short_term_commands(Node& node, CmdManager& manager)
    {
        auto cmdsMap = manager.getCurrentTypeToStr();
        std::unordered_map<Cmd, std::string> cmd_map;
        std::string prefix_cmd = "read_aiger " + node.aig_addr.string() + "; strash";

        fs::path st_save_path = node.aig_addr.parent_path();
        
        bool have_his_st_file = false;
        for (const auto& entry : fs::directory_iterator(st_save_path)) {
            fs::file_status fstatus = entry.status();  
            if (fs::is_regular_file(fstatus) && entry.path().filename().string().find("cb_st_") == 0) {
                have_his_st_file = true; 
            }
        }
        
        for (const auto& c : cmdsMap) {
            Cmd arm = c.first;
            
            std::string cmd = "";
            if (have_his_st_file){
                cmd = fmt::format("read_aiger {}/cb_st_{}.aiger; strash; ", st_save_path.string(), static_cast<int>(arm));
            } else{
                cmd = fmt::format("{}; {}; strash; write_aiger {}/cb_st_{}.aiger; strash; ", prefix_cmd, c.second, st_save_path.string(), static_cast<int>(arm));   
            }

            cmd += (target >= 2) ? ((exp_num == 1) ? fpga_cmd_flowtune : fpga_cmd) : "";
            cmd += "print_stats; ";
            cmd_map[arm] = cmd; 
        }
        
        return cmd_map;
    }

    std::string constrain_random_opts(CmdManager& manager, int search_depth) {
        int numArms = manager.getCurrentMapSize();
        auto cmdsMap = manager.getCurrentTypeToStr();

        std::string random_seq = "";
        std::random_device rd;
        std::mt19937 gen(rd());

        for (int i = 0; i < search_depth; i++) {
            std::uniform_int_distribution<> distrib(0, numArms - 1);

            auto iter = cmdsMap.begin();
            std::advance(iter, distrib(gen));  

            random_seq += iter->second + "; ";
        }
        return random_seq;
    }

    std::unordered_map<Cmd, std::vector<std::string>> constr_long_term_commands(Node& node, CmdManager& manager, int search_times, const std::vector<int>& search_depths)
    {
        int numArms = manager.getCurrentMapSize();
        auto cmdsMap = manager.getCurrentTypeToStr();

        std::unordered_map<Cmd, std::vector<std::string>> cmd_map;
        std::string prefix_cmd = "read_aiger " + node.aig_addr.string() + "; strash";
        fs::path lt_save_path = node.aig_addr.parent_path();
 
        for (const auto& c : cmdsMap) {
            Cmd arm = c.first;

            for (int j = 0; j < search_depths.size(); j++){
                int depth = search_depths[j];

                for (int i = 0; i < search_times; i++){

                    std::string cmd = fmt::format("{}; {}; {}strash; write_aiger {}/cb_lt_{}_{}.aiger; strash; ", prefix_cmd, c.second, constrain_random_opts(manager, depth-1), lt_save_path.string(), static_cast<int>(arm), j*search_times+i); 

                    cmd += (target >= 2) ? ((exp_num == 1) ? fpga_cmd_flowtune : fpga_cmd) : "";
                    cmd += "print_stats; ";
                    cmd_map[arm].push_back(cmd);
                }
            }
        }
        return cmd_map;
    }

    Node next_rand(CmdManager& manager, std::vector<Cmd> types, int l = 1, int r = 3, int suffix = 0) const {
        std::vector<Cmd> types_noBL;
        for (auto t : types){
            if (t != Cmd::B)  
                types_noBL.push_back(t);
        }

        fs::path res_next_aig = aig_addr.parent_path() / "children"/("res" + std::to_string(suffix) + ".aiger");

        int len = rand_int(l,r);
        std::vector<Oper> res_seq = seq;
        int node_length = length;
        std::string cmd = "read_aiger " + aig_addr.string() + "; strash; ";
        
        for (int i = 0; i < len; ++i){
            Oper oper;
            if (res_seq.size() && (res_seq.back().type == Cmd::B)){
                oper = rand_oper(types_noBL);
            } else 
                oper = rand_oper(types);
            
            res_seq.emplace_back(oper);
            node_length++;
            cmd += oper2str(manager, oper)+"; ";
        }    

        Node res_node; 
        res_node.seq = res_seq;
        res_node.length = node_length;
        res_node.aig_addr = res_next_aig;

        if (target == 0 || target == 1){
            cmd += fmt::format("write_aiger {}; strash; &get; &ps; ", res_next_aig.string());
            aig_res ress = res_aig(remove_ansi_sequences(executeAbcCommand(cmd)));
            res_node.aig_node = ress.aig_node;
            res_node.aig_lev = ress.levs;
        } else{
            cmd += fmt::format("write_aiger {}; strash; {}print_stats; ", res_next_aig.string(), exp_num == 1 ? fpga_cmd_flowtune : fpga_cmd);
            map_res ress = res_lut(executeAbcCommand(cmd));
            res_node.lut_area = ress.lut_nd;
            res_node.lut_lev = ress.lut_lev;
        }
        
        return res_node;
    }

    Node(const Node& other)
        : aig_addr(other.aig_addr),
        seq(other.seq),
        length(other.length),
        aig_node(other.aig_node),
        aig_lev(other.aig_lev),
        lut_area(other.lut_area),
        lut_lev(other.lut_lev)
    {}

    Node(){}
};

struct less_node {
    /**
     * Compare two nodes based on their optimization targets.
     * For AIG targets (0,1): compares AIG nodes and levels
     * For LUT targets (2,3): compares LUT area and levels
     */
    less_node() = default;

    bool operator()(const Node &n0, const Node &n1) const {
        // Helper functions to calculate costs
        auto calc_aig_cost = [](const Node &n) -> long {
            return target == 0 ? n.aig_node : n.aig_node * n.aig_lev;
        };
        
        auto calc_lut_cost = [](const Node &n) -> long {
            return target == 2 ? n.lut_area : n.lut_area * n.lut_lev;
        };

        // Compare based on target type
        if (target < 2) {  // AIG optimization
            const long cost0 = calc_aig_cost(n0);
            const long cost1 = calc_aig_cost(n1);
            
            if (cost0 != cost1) return cost0 < cost1;
            if (n0.aig_node != n1.aig_node) return n0.aig_node < n1.aig_node;
        } else {  // LUT optimization
            const long cost0 = calc_lut_cost(n0);
            const long cost1 = calc_lut_cost(n1);
            
            if (cost0 != cost1) return cost0 < cost1;
            if (n0.lut_lev != n1.lut_lev) return n0.lut_lev < n1.lut_lev;
        }

        // If all other metrics are equal, compare sequence length
        return n0.seq.size() < n1.seq.size();
    }
};

        

struct eq_node {
    eq_node() {}

    bool operator()(const Node &n0, const Node &n1) const {
        if (target == 0 || target == 1) {
            return n0.aig_node == n1.aig_node && n0.aig_lev == n1.aig_lev;
        } else {
            return n0.lut_area == n1.lut_area && n0.lut_lev == n1.lut_lev;
        }
    }
};

Node get_ori_node(const fs::path& i_path){
    Node ori_node; 

    std::string extension = i_path.extension().string();
    std::string read_cmd = "";
    if (extension == ".aig" || extension == ".aiger") {
        read_cmd = "read_aiger";
    } else if (extension == ".v") {
        read_cmd = "read_verilog";
    } else if (extension == ".blif") {
        read_cmd = "read_blif";
    } else {
        throw std::runtime_error("Unsupported file extension.");
    }

    // Get the original aig nodes/levels
    fs::create_directory(tmp_dir_path);
    fs::path dest_path = tmp_dir_path / "origin.aiger";

    std::string cmd = fmt::format("{} {}; strash; print_stats; write_aiger {};", read_cmd, i_path.string(), dest_path.string());
    aig_res ress_orig = res_aig(remove_ansi_sequences(executeAbcCommand(cmd)));
    ori_node.aig_addr = dest_path;
    ori_node.aig_node = ress_orig.aig_node;
    ori_node.aig_lev = ress_orig.levs;
    set_time_limit_by_aig_nodes(ress_orig.aig_node);

    std::cout<<"original :"<<std::endl;
    // Get the fpga mapping results
    if (target >= 2){
        if (exp_num == 1 ){
            cmd = fmt::format("{} {}; strash; {}print_stats; ", read_cmd, i_path.string(), fpga_cmd_flowtune);
        } else{
            cmd = fmt::format("{} {}; strash; {}print_stats; ", read_cmd, i_path.string(), fpga_cmd);
        }

        map_res ress_orig = res_lut(remove_ansi_sequences(executeAbcCommand(cmd)));
        ori_node.lut_area = ress_orig.lut_nd;
        ori_node.lut_lev = ress_orig.lut_lev;
        set_time_limit_by_aig_nodes(ress_orig.aig_node);

        std::cout<<"lut  area: "<<ori_node.lut_area<<"    lut level: "<<ori_node.lut_lev << std::endl;
    } else{
        std::cout<<"aig  node: "<<ori_node.aig_node<<"    aig level: "<<ori_node.aig_lev << std::endl;
    }

    return ori_node;
}

Cmd stringToCmd(const std::string& cmd) {
    if (cmd == "rewrite") return Cmd::RW;
    if (cmd == "rewrite -z") return Cmd::RWZ;
    if (cmd == "balance") return Cmd::B;
    if (cmd == "refactor") return Cmd::RF;
    if (cmd == "refactor -z") return Cmd::RFZ;
    if (cmd == "resub") return Cmd::RESUB;
    if (cmd == "resub -z") return Cmd::RESUBZ;
    // Other commands
    throw std::invalid_argument("Unsupported or unknown command: " + cmd);
}

Oper str2oper(const std::string& soper){
    Oper oper;
    std::vector<std::string> tokens = split(soper,' ');
    oper.type = stringToCmd(tokens[0]);

    switch (oper.type) {
        case Cmd::RW:
            oper.f0 = false;
            for (size_t i = 1; i < tokens.size(); ++i) {
                if (tokens[i] == "-z") {
                    oper.f0 = true;
                    oper.type = Cmd::RWZ;
                }
            }
            break;
        case Cmd::RF:
            oper.f0 = false;
            for (size_t i = 1; i < tokens.size(); ++i) {
                if (tokens[i] == "-z") {
                    oper.f0 = true;
                    oper.type = Cmd::RFZ;
                }
            }
            break;
        case Cmd::RESUB:
            oper.f0 = false;
            for (size_t i = 1; i < tokens.size(); ++i) {
                if (tokens[i] == "-z") {
                    oper.f0 = true;
                    oper.type = Cmd::RESUBZ;
                }
            }
            break;
        case Cmd::B:
            oper.f0=false;
            break;
        // Other commands
        default:
            break;
    }
    return oper;
}

// Use for sequence output
std::string oper2str(CmdManager& manager, const Oper &oper) {

    std::string res = "";
    res += manager.getFullCommand(oper.type);

    switch (oper.type) {
        case Cmd::SPECIAL_AIG:
            res += add_aig_opt_cmd;
            break;
        case Cmd::SPECIAL_LUT:
            res += add_fpga_opt_cmd;
            break;
        default:
            break;
    }
    return res;
}

void write_ans(const Node &best_node, const std::string& path_out, CmdManager& manager){
    std::ofstream fout(path_out);
    fout << best_node.get_seq_str(manager) << std::endl;
    fout.close();
}

// Random parameters for operations
Oper rand_oper(std::vector<Cmd> types){
    Cmd type = types[rand()%types.size()];
    Oper oper(type);
    return oper;
}

// Random parameters for operations
// Oper rand_oper(std::vector<Cmd> types){
//     Cmd type = types[rand()%types.size()];
//     Oper oper(type);

//     switch (type) {
//         case Cmd::RW:
//             oper.f0 = rand()%2;
//             break;
//         case Cmd::RWZ:
//             oper.f0 = 1;
//             oper.f1 = rand()%2;      
//             break; 
//         case Cmd::B:
//             oper.f0 = rand()%2;
//             break;
//         case Cmd::RF:
//             oper.f0 = rand()%2;
//             oper.f1 = rand()%2;
//             break;
//         case Cmd::RFZ:
//             oper.f0 = rand()%2;
//             oper.f1 = rand()%2;
//             break;
//         case Cmd::RESUB:
//             oper.arg0 = rand_int(6,12);
//             break;
//         case Cmd::RESUBZ:
//             oper.arg0 = rand_int(6,12);
//             break;
//         default:
//             // Optionally handle unexpected types or do nothing
//             break;
//     }
//     return oper;
// }

