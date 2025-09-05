#ifndef _EVOLUTION_
#define _EVOLUTION_

#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <climits>
#include <filesystem>
#include <omp.h>
#include <algorithm>
#include <fmt/core.h>
#include "utils.hpp"
#include "cmd.hpp"
#include "node.hpp"

namespace fs = std::filesystem;

namespace Evolution{

// All function declarations 
void report_res(const Node& node, CmdManager& manager);
Node evo_search(CmdManager& manager, const Node& stage_node, const std::string& path_out, uint32_t limit_population, uint32_t select_rate);

void report_res(const Node& node, CmdManager& manager) {
    if (target == 0){
        std::cout << "    aig  node: " << node.aig_node << "    aig level: " << node.aig_lev;
    } else{
        std::cout << "    lut  area: " << node.lut_area << "    lut level: " << node.lut_lev;
    }
    std::cout<<"    len seq: "<< node.seq.size()<<"    seq: "<< node.get_seq_str(manager)<<std::endl;
}

// Seqfinder: search the synthesis flow
Node evo_search(CmdManager& manager, const Node& stage_node, const std::string& path_out, uint32_t limit_population, uint32_t select_rate){

    // Set random seed for reproducibility
    unsigned seed = 12345;
    std::srand(seed);

    int seq_length_evo = 0;
    if (seq_all != 0){
        if (exp_num == 1) {
            seq_length_evo = seq_all - 14;
        } else {
            seq_length_evo = seq_all - static_cast<int>(seq_all / 4);
        }
        time_limit = INT_MAX;
        std::cout << "Evolution tuning using seq length limit to: " << seq_length_evo << std::endl;
    } else{
        seq_length_evo = INT_MAX;
        std::cout << "Evolution tuning using time limit (s): " << time_limit << std::endl;
    }

    tic_toc t_c;
    
    Node best_node = stage_node;
    std::cout<<"    original :"<<std::endl;
    report_res(best_node, manager);

    // Initial population
    int total_size = static_cast<int>(limit_population * select_rate);
    std::vector<Node> population(total_size);
    fs::create_directories(best_node.aig_addr.parent_path()/ "children");

    #pragma omp parallel for
    for (int i = 0; i < total_size; ++i){
        Node result;
        result = best_node.next_rand(manager, {Cmd::RW, Cmd::RWZ, Cmd::B, Cmd::RF, Cmd::RFZ, Cmd::RESUB, Cmd::RESUBZ}, 1, 3, i);
        
        #pragma omp critical
        {
            population[i] = result;
        }
    }

    // Choice command settings
    int choice_cnt = 0;
    int choice_time = time_limit / 2;
    int choice_length = seq_length_evo / 2;
    bool choice_flag = false;
    bool choice_compute_flag = false;

    int converge_cnt = 0;
    for(int idx_gen = 0; idx_gen < population.size() && t_c.toc() <= time_limit && best_node.length <= seq_length_evo; ++idx_gen){
        std::cout << "---------------------------------------------------";
        std::cout << "\n    generation: " << idx_gen << "    time(s): " << t_c.toc() << "    population size: " << population.size() << std::endl;

        //----------selection----------
        sort(population.begin(), population.end(), less_node()); 
        if (population.size() > limit_population){
            population.resize(limit_population);
        }

        if (!choice_flag){
            if (
                ((best_node.length < choice_length || t_c.toc() <= choice_time) && (std::rand() % 3)) ||
                ((best_node.length >= choice_length || t_c.toc() > choice_time) && !choice_compute_flag)
            ) {
                choice_cnt++;
                for(size_t j = 0; j < population.size(); ++j) {
                    auto &node = population[j];
                    node.seq.push_back(rand_oper({Cmd::CHOICE_S}));
                
                    if (choice_cnt == 3){
                        node.seq.push_back(rand_oper({Cmd::CHOICE_C}));
                        choice_flag = true;
                    }
                }

                if (choice_cnt == 3){
                    choice_cnt = 0;
                    choice_compute_flag = true;
                } else{
                    choice_compute_flag = false;
                }
            }
        }
    
        if (less_node()(population.front(), best_node)) {
            best_node = population.front();
            write_ans(best_node, path_out, manager);
            report_res(best_node, manager);
            converge_cnt = 0;   
        } else{
            converge_cnt++;
        }

        if ((converge_cnt > 5) && (!choice_compute_flag)){
            break;
        }
        //------------------------------------
        
        //---------next_generation------------
        std::vector<Node> next_population;
        static int global_cnt = 0; 
        #pragma omp parallel  
        {
            std::vector<Node> local_next_population;  

            #pragma omp for
            for(size_t j = 0; j < population.size(); ++j) {
                auto node = population[j];  // Remove const to allow modification
                
                // Add suffix sequence for nodes with length multiple of 14
                if (exp_num == 1 && node.length % 14 == 0){
                    if (target == 0 || target == 1){
                        std::string cmd = fmt::format("read_aiger {}; {}; write_aiger {}; print_stats;", node.aig_addr.string(), add_aig_opt_cmd, node.aig_addr.string());
                        aig_res ress = res_aig(remove_ansi_sequences(executeAbcCommand(cmd)));
                        node.aig_node = ress.aig_node;
                        node.aig_lev = ress.levs;
                        node.seq.push_back(Oper(Cmd::SPECIAL_AIG));
                    } else {
                        std::string cmd = fmt::format("read_aiger {}; {}; write_aiger {}; if -K 6; print_stats;", node.aig_addr.string(), add_fpga_opt_cmd, node.aig_addr.string());
                        map_res ress = res_lut(remove_ansi_sequences(executeAbcCommand(cmd)));
                        node.lut_area = ress.lut_nd;
                        node.lut_lev = ress.lut_lev;
                        node.seq.push_back(Oper(Cmd::SPECIAL_LUT));
                    }
                }
                
                fs::create_directories(node.aig_addr.parent_path()/ "children");
                bool compressed = false;

                for(auto i = 0u; i < select_rate; ++i) {
                    if (i == 3 && !compressed){
                        break;
                    }
                        
                    int local_cnt;
                    #pragma omp atomic capture
                    local_cnt = global_cnt++;  
                    Node next_node = node.next_rand(manager, {Cmd::RW, Cmd::RWZ, Cmd::B, Cmd::RF, Cmd::RFZ, Cmd::RESUB, Cmd::RESUBZ}, 1, 2, local_cnt);
                    if (!eq_node()(next_node, node)) {
                        local_next_population.emplace_back(next_node);
                        compressed = true;
                    }
                }
            }

            #pragma omp critical
            {
                next_population.insert(next_population.end(), local_next_population.begin(), local_next_population.end());
            }
        }

        population = std::move(next_population);
        //------------------------------------
    }

    // if (exp_num == 1){
    //     if (target <= 1){
    //         if (best_node.seq.empty() || best_node.seq.back().type != Cmd::SPECIAL_AIG) {
    //             std::string cmd = fmt::format("read_aiger {}; {}; write_aiger {}; print_stats;", best_node.aig_addr.string(), add_aig_opt_cmd, best_node.aig_addr.string());
    //             aig_res ress = res_aig(remove_ansi_sequences(executeAbcCommand(cmd)));
    //             best_node.aig_node = ress.aig_node;
    //             best_node.aig_lev = ress.levs;
    //             best_node.seq.push_back(Oper(Cmd::SPECIAL_AIG));
    //         }
    //     } else {
    //         if (best_node.seq.empty() || best_node.seq.back().type != Cmd::SPECIAL_LUT) {
    //             std::string cmd = fmt::format("read_aiger {}; {}; write_aiger {}; if -K 6; print_stats;", best_node.aig_addr.string(), add_fpga_opt_cmd, best_node.aig_addr.string());
    //             map_res ress = res_lut(remove_ansi_sequences(executeAbcCommand(cmd)));
    //             best_node.lut_area = ress.lut_nd;
    //             best_node.lut_lev = ress.lut_lev;
    //             best_node.seq.push_back(Oper(Cmd::SPECIAL_LUT));
    //         }
    //     }
    // }
    write_ans(best_node, path_out, manager);

    return best_node;
}

}

#endif