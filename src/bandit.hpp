#pragma once

#include <vector>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <unordered_map>
#include <cmath>
#include <omp.h>
#include <algorithm>
#include <random>
#include <fmt/core.h>
#include <numeric>
#include <Eigen/Dense>
#include "cmd.hpp"
#include "utils.hpp"
#include "features.hpp"
#include "node.hpp"

namespace fs = std::filesystem;
namespace Bandit{

// Default settings
tic_toc t_c;
int search_times = 3;
int cb_iter = 10;
double decayRate = 0.4;
int seq_bandit_length = 0; // Default value is 0, which would be set in main.cpp

// Stage info for return-back
struct StageInfo {
    Node node;          // Result node for this stage
    Cmd cmd;            // Command used in this stage
    bool needs_retry;   // Whether this stage needs to be retried
    
    StageInfo() : needs_retry(false) {}  // Default constructor
    void clear() {
        node = Node();
        cmd = Cmd::NONE;
        needs_retry = false;
    }
};

// Function definitions
Node get_ori_node(const fs::path& i_path);
Cmd Agent(Node& node, CmdManager& manager, int stage, HashTable& valuetable, int* curr_cost);
std::unordered_map<Cmd, std::vector<float>> featureScale(std::unordered_map<Cmd, std::vector<float>> data, float scale_ratio);
std::unordered_map<Cmd, float> scaleData(std::unordered_map<Cmd, Node> node_st_res, float scale_ratio);

// Bandit algorithm definition
class CBandit {
public:
    CBandit(CmdManager& manager, std::vector<Cmd>& ar, unsigned int dim, const std::vector<double>& weights)
        : manager(manager), arms(ar), d(dim), rho(0.1), w_vector(Eigen::VectorXd::Map(weights.data(), weights.size())), w_sum(w_vector.sum()), alpha_val(1.0 + sqrt(log(2.0 / rho) / 2.0)) {
        for (auto arm : arms) {
            a_matrix.push_back(Eigen::MatrixXd::Identity(d, d));
            b_vector.push_back(Eigen::VectorXd::Zero(d));
            a_matrix_inverse.push_back(Eigen::MatrixXd::Identity(d, d)); // Initialize with identity matrix
        }
        p_arm.resize(arms.size(), INIT_SCORE);

    }

    void set_rho(double r) {
        rho = r;
        alpha_val = 1.0 + sqrt(log(2.0 / rho) / 2.0);
    }

    void score_arm(Eigen::VectorXd& x, int index, int chosen_times) {
        double p = payoff(x, index, chosen_times);
        p_arm[index] = p;
    }

    Cmd best_arm() {
        double max_value = *std::max_element(p_arm.begin(), p_arm.end());

        // Find the index of the maximum value
        std::vector<int> max_indices;
        for (int i = 0; i < p_arm.size(); ++i) {
            if (p_arm[i] == max_value) {
                max_indices.push_back(i);
            }
        }

        // If there is only one maximum, return that arm
        if (max_indices.size() == 1) {
            return arms[max_indices[0]];
        }

        // If there are multiple maxima, choose a random arm
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, max_indices.size() - 1);
        int random_index = max_indices[dis(gen)];

        return arms[random_index];
    }

    void print_score() {
        std::cout << "Each arm and its corresponding score: " << std::endl;

        for (size_t i = 0; i < arms.size(); ++i) {
            std::cout << manager.getFullCommand(arms[i]) << ": " << p_arm[i] << std::endl;
        }

    }

    void update(int i, Eigen::VectorXd& x, double reward) {
        // double weight = w_vector(i) / w_sum;
        a_matrix[i] += x * x.transpose();
        b_vector[i] += reward * x;

        // Update cached inverse matrix
        a_matrix_inverse[i] = a_matrix[i].inverse();
    }

private:
    double softmax(double value, double sum, double temperature) {
        return exp(value / temperature) / sum;
    }

    double payoff(Eigen::VectorXd& context, int index, int chosen_times) {
        Eigen::VectorXd weighted_context = context.array() * w_vector.array();
        Eigen::VectorXd o = a_matrix_inverse[index] * b_vector[index];
        Eigen::VectorXd t = weighted_context.transpose() * a_matrix_inverse[index];

        // Multiply each element of the context vector by its corresponding weight
        alpha_val = 1.0 + sqrt(log(2.0 / rho) / chosen_times);

        double p = o.transpose().dot(weighted_context) + alpha_val * sqrt(t.dot(weighted_context));

        return p;
    }

    std::vector<Cmd> arms;
    unsigned int d;
    double rho;
    std::vector<Eigen::MatrixXd> a_matrix;
    std::vector<Eigen::MatrixXd> a_matrix_inverse;
    std::vector<Eigen::VectorXd> b_vector;
    std::vector<double> p_arm;
    Eigen::VectorXd w_vector;
    double w_sum;
    double alpha_val;
    CmdManager& manager;

    static const double INIT_SCORE;
};

const double CBandit::INIT_SCORE = 0.0;

std::unordered_map<Cmd, std::vector<float>> featureScale(std::unordered_map<Cmd, std::vector<float>> data, float scale_ratio) {
    std::unordered_map<Cmd, std::vector<float>> scaled_data;

    // Get the sum of all the data
    size_t expectedLength = data.begin()->second.size();
    std::vector<float> sums(expectedLength, 0.0f);
    std::vector<std::vector<float>> tmp_values(expectedLength, std::vector<float>());
    std::vector<bool> allEqual(expectedLength, true);
    std::vector<float> firstValues(expectedLength);

    bool isFirst = true;
    for (const auto& pair : data) {
        const std::vector<float>& vec = pair.second;
        if (isFirst) {
            std::copy(vec.begin(), vec.end(), firstValues.begin());
            isFirst = false;
        }

        for (size_t i = 0; i < vec.size(); ++i) {
            sums[i] += vec[i];
            tmp_values[i].push_back(vec[i]);
            if (vec[i] != firstValues[i]) {
                allEqual[i] = false;
            }
        }
    }

    for (size_t i = 0; i < expectedLength; ++i) {
        if (allEqual[i]) {
            for (const auto& pair : data) {
                scaled_data[pair.first].push_back(1.0f);
            }
        } else{
            float average = sums[i] / data.size();
            std::vector<float> tmp_value = tmp_values[i];

            float maxVal = std::numeric_limits<float>::lowest();  
            float minVal = std::numeric_limits<float>::max();

            for (int k = 0; k < tmp_value.size(); k++){
                float scaledValue = pow(average / tmp_value[k], scale_ratio);
                tmp_value[k] = scaledValue;

                // Find the max and min value of the vector
                if (scaledValue > maxVal) {
                    maxVal = scaledValue;
                }
                if (scaledValue < minVal) {
                    minVal = scaledValue;
                }
            }

            float range = maxVal - minVal;
            for (int k = 0; k < tmp_value.size(); k++) {
                tmp_value[k] =  (tmp_value[k] - minVal) / range * 2.0f + 0.001f;
            }

            int k = 0;
            for (const auto& pair : data) {
                scaled_data[pair.first].push_back(tmp_value[k]);
                k++;
            }
        }
    }

    return scaled_data;
}

std::unordered_map<Cmd, float> scaleData(std::unordered_map<Cmd, Node> node_st_res, float scale_ratio) {
    std::unordered_map<Cmd, float> scaled_st;
    std::unordered_map<Cmd, int> tmp_cost;
    
    // Get the sum of all the data
    int sum = 0;
    bool allSame = true;

    auto it = node_st_res.begin();
    const Cmd& first_cmd = it->first;
    const Node& first_node = it->second;
    int first_cost = first_node.calculate_cost();
       
    sum += first_cost;
    tmp_cost[first_cmd] = first_cost;

    ++it; 
    for (; it != node_st_res.end(); ++it) {
        const Cmd& cmd = it->first;
        const Node& tmp_node = it->second;
        int cost = tmp_node.calculate_cost();

        if (allSame && cost != first_cost) 
            allSame = false;

        sum += cost;
        tmp_cost[cmd] = cost;
    }
   
    // Scale the data
    if (allSame) {
        // std::cout << "Data fully same!"<<std::endl;
        for (const auto& pair : node_st_res) {
            scaled_st[pair.first] = 1.0f;
        }

    } else {
        float average = static_cast<float>(sum / node_st_res.size());

        // Scale each value based on the average
        for (const auto& pair : tmp_cost) {
            float scaledValue = average / pair.second;
            scaledValue = pow(scaledValue, scale_ratio);
            scaled_st[pair.first] = scaledValue;
        }

        auto extremalVal = findExtremalForId(scaled_st);
        float range = extremalVal.first - extremalVal.second;
        
         for (const auto& pair : scaled_st) {
            float extendDiffVal = (pair.second - extremalVal.second) / range * 2.0f;
            scaled_st[pair.first] = extendDiffVal + 0.001f;
        }

    }

    return scaled_st;
}

void updateContext(Node& node, CmdManager& manager, const std::vector<int>& search_depths, int stage,  
                    LongTermStatsMap<float>& long_term_avg_map, LongTermStatsMap<int>& long_term_min_map, 
                    std::unordered_map<Cmd, Eigen::VectorXd>& context_info, int d, HashTable& valuetable, int initial_flag)
{
    // Default settings
    int numArms = manager.getCurrentMapSize();
    auto cmdsMap = manager.getCurrentTypeToStr();
    if (initial_flag == 0)
        search_times = 1;

    std::unordered_map<Cmd, std::vector<Node>> node_lt_res = node.get_long_term_reward(node, manager, search_times, search_depths);

    // Calculate long-term average and minimum for each arm
    std::unordered_map<Cmd, std::vector<int>> min_val;
    std::unordered_map<Cmd, std::vector<int>> sum_val;
    
    for (const auto& c : cmdsMap) {
        Cmd arm = c.first;

        min_val[arm] = std::vector<int>(search_depths.size(), INT_MAX);
        sum_val[arm] = std::vector<int>(search_depths.size(), 0.0);

        for (int i = 0; i < search_depths.size(); i++) {
            int search_depth = search_depths[i];
            for (int j = 0; j < search_times; ++j) {
                int index = i*search_times + j;
                Node tmp_node = node_lt_res[arm][index];
                int val = tmp_node.calculate_cost();               

                // Update the depth informatio of each arm
                sum_val[arm][i] += val;
                if (val < min_val[arm][i]) {
                    min_val[arm][i] = val;
                }
            }

            if (initial_flag == 0){ 
                float his_val_avg = long_term_avg_map.get(search_depth, arm);
                long_term_avg_map.update(search_depth, arm, (his_val_avg+sum_val[arm][i])/(search_times+1));
                int his_val_min = long_term_min_map.get(search_depth, arm);
                long_term_min_map.update(search_depth, arm, std::min(his_val_min, min_val[arm][i]));
            }else{
                // Corresponding to the initial situation
                long_term_avg_map.update(search_depth, arm, sum_val[arm][i] / search_times);
                long_term_min_map.update(search_depth, arm, min_val[arm][i]);
            }             
        }
    }

    for (int i = 0; i < search_depths.size(); ++i) {
        int search_depth = search_depths[i];

        // Update the hashtable and get the min value of all arms
        int min_value = INT_MAX;
        min_value = std::min(min_value, long_term_min_map.findMinForId(search_depth));
        valuetable[stage][search_depth+stage-1] = min_value;

        // Scale data of the long_term_min_map
        long_term_avg_map.scaleData(search_depth, 8.0);
        long_term_min_map.scaleData(search_depth, 8.0);
    }

    // Update context_info
    for (const auto& c : cmdsMap) {
        Cmd arm = c.first;

        int offset = d; 
        for (int i = 0; i < search_depths.size(); ++i) {
            int search_depth = search_depths[i];

            context_info[arm][offset] = long_term_avg_map.getScaledData(search_depth, arm);
            context_info[arm][offset + 1] = long_term_min_map.getScaledData(search_depth, arm);
            offset += 2;
        }
    }
}

int modifyIterations(int cb_iter, int stage, int numArms) {
    int minIterations = numArms;
    int rlIterations = cb_iter - static_cast<int>(decayRate * stage);
    // rlIterations = cb_iter;
    rlIterations = std::max(rlIterations, minIterations);
    return rlIterations;
}

Cmd Agent(Node& node, CmdManager& manager, int stage, HashTable& valuetable, int* curr_cost)
{
    int numArms = manager.getCurrentMapSize();
    auto cmdsMap = manager.getCurrentTypeToStr();
    int iter = modifyIterations(cb_iter, stage, numArms);

    //////////////////Static features extractions////////////////
    // Store the short term reward
    std::unordered_map<Cmd, Node> node_st_res = node.get_short_term_reward(node, manager);
    Node st_best_node;
    Cmd st_ori_cmd;
    if (!node_st_res.empty()) {
        auto it = node_st_res.begin(); 
        st_ori_cmd = it->first;
        st_best_node = it->second;
    }

    for (const auto& pair : node_st_res) {
        if (pair.first == st_ori_cmd)
            continue;
        if (less_node()(pair.second, st_best_node)){
            st_best_node = pair.second;
        }  
    }

    valuetable[stage][stage] = st_best_node.calculate_cost();

    int d = 2; // feature dimension
    std::unordered_map<Cmd, std::vector<float>> static_features;  
    std::vector<std::pair<Cmd, std::string>> cmdsVec(cmdsMap.begin(), cmdsMap.end());
    std::unordered_map<Cmd, std::vector<float>> local_features;
    int use_ccirc = exp_num == 1 ? 0 : 1; // 0: use yosys, 1: use ccirc 

    #pragma omp parallel for
    for (size_t i = 0; i < cmdsVec.size(); ++i) {
        Cmd arm = cmdsVec[i].first;

        std::vector<float> extracted_data;
        if (use_static_extraction) {  // Use static analysis tools (ccirc/yosys)
            extracted_data = extractInputData(node_st_res[arm], stage, arm, use_ccirc);
        } 

        Node tmp_node = node_st_res.at(arm);
        if (target == 0 || target == 1){
            extracted_data.push_back(tmp_node.aig_node);
            extracted_data.push_back(tmp_node.aig_lev);
        } else {
            extracted_data.push_back(tmp_node.lut_area);
            extracted_data.push_back(tmp_node.lut_lev);
        }

        #pragma omp critical
        {
            local_features[arm] = extracted_data; 
        }
    }

    if (use_static_extraction && use_ccirc) {
        std::string rm_tmp_blif = fmt::format("rm netlist_*");
        system(rm_tmp_blif.c_str()); 
    }
    
    static_features = featureScale(local_features, 8.0);

    //////////////////Initialize the contextual information////////////////
    std::vector<int> search_depth_list;  
    std::vector<double> weights;

    if (use_static_extraction) {
        if (use_ccirc) {
            weights = {0.083,0.083,0.014,0.089,0.121,0.121,0.121,0.1,0.1,0.1,0.1,0.1,0.1,0.021,0.021,0.021,0.021,0.021,0.008,0.008, 0.1, 0.1};
            d += 22;
        } else {
            weights = {0.083,0.083,0.083};
            d += 3;
        }
    }
    if (stage < static_cast<int>(seq_bandit_length / 3)) {
        search_depth_list = {2,3,5};
        weights.insert(weights.end(), {0.4,0.2,0.5,0.55,0.8,0.85,1.3,1.35});
    } else if (stage < static_cast<int>(seq_bandit_length * 2 / 3)) {
        search_depth_list = {2,3,4};
        weights.insert(weights.end(), {0.4,0.2,0.5,0.55,0.75,0.8,0.95,1.0});
    } else {
        search_depth_list = {2,3};
        weights.insert(weights.end(), {0.4,0.2,0.5,0.55,0.75,0.8});
    } 

    // Check if weight dimension matches feature dimension
    assert((weights.size() == d + 2*search_depth_list.size()) && "Error: weights.size() does not match feature dimension (d + 2*search_depth_list.size())");

    // Contextual bandit initialization
    std::vector<Cmd> arms = manager.getArms();
    std::unordered_map<Cmd, int> armChosenTimes;
    LongTermStatsMap<float> long_term_avg_map; // average results of all arms for each depth
    LongTermStatsMap<int> long_term_min_map; // min results of all arms for each depth
    std::unordered_map<Cmd, Eigen::VectorXd> context_info;
    
    for (const auto& c : cmdsMap) {
        Cmd arm = c.first;

        Eigen::VectorXd context(d + 2*search_depth_list.size());
        for (int j = 0; j < d; ++j){
            context(j) = static_features[arm][j];
        }
        
        int offset = d;
        for (int depth_index = 0; depth_index < search_depth_list.size(); ++depth_index){
            context(offset) = 0;
            context(offset + 1) = 0;
            offset += 2;
        }
        context_info[arm] = context;
        armChosenTimes[arm] = 1;
    }    
    
    CBandit* C = new CBandit(manager, arms, d + 2*search_depth_list.size(), weights); 
    C->set_rho(0.2);  

    updateContext(node, manager, search_depth_list, stage, long_term_avg_map, long_term_min_map, context_info, d, valuetable, 1);
    
    //////////////////Add dynamic features during iterations////////////////
    std::cout<<"Start iteration!"<<std::endl;
    bool break_loop = false;
    Cmd chosen;      // Current chosen command
    Cmd prev_chosen; // Last chosen command
    int consecutive_same_chosen = 0; 
    for (int t = 1; t <= iter; t++){
        // Update the context information (long term results)
        updateContext(node, manager, search_depth_list, stage, long_term_avg_map, long_term_min_map, context_info, d,valuetable, 0);

        for (const auto& c : cmdsMap) {
            Cmd arm = c.first;
            C->score_arm(context_info[arm], manager.getCmdIndex(arm), armChosenTimes[arm]);
        } 

        C->print_score();
        chosen = C->best_arm(); 
        armChosenTimes[chosen] ++;
        std::cout << fmt::format("iter: {}; chosen: {};", t, cmdsMap[chosen]) << std::endl;

        // Check if three consecutive choices are the same
        if (t != 1){
            if (chosen == prev_chosen) {
                consecutive_same_chosen++;
            } else {
                consecutive_same_chosen = 0;
            }

            std::cout << "Consecutive Choices: "<<consecutive_same_chosen<<std::endl;
            prev_chosen = chosen;
        }
        
        if (consecutive_same_chosen >= 2) {
            std::cout << "Consecutive same chosen, exiting loop." << std::endl;
            break_loop = true;
        }

        if (break_loop) {
            break;
        }
        
        std::unordered_map<Cmd, float> scaled_st = scaleData(node_st_res, static_cast<float>(t));

        float r = scaled_st[chosen];
        C->update(manager.getCmdIndex(chosen), context_info[chosen], r);
        std::cout << "--------------------------------------" << std::endl;
    }

    Node curr_cost_node = node_st_res[chosen];
    *curr_cost = curr_cost_node.calculate_cost();

    node.seq.push_back(Oper(chosen));
    if (target == 0 || target == 1){
        node.aig_node = curr_cost_node.aig_node;
        node.aig_lev = curr_cost_node.aig_lev;
    } else {
        node.lut_area = curr_cost_node.lut_area;
        node.lut_lev = curr_cost_node.lut_lev;
    }

    fs::path next_aig_addr = node.aig_addr.parent_path() / "children" / ("res" + std::to_string(stage) + ".aiger");
    fs::create_directories(next_aig_addr.parent_path());
    std::string cmd = fmt::format("read_aiger {}; {}; strash; write_aiger {};", node.aig_addr.string(), manager.getFullCommand(chosen), next_aig_addr.string());
    std::string output = executeAbcCommand(cmd);
    node.aig_addr = next_aig_addr;
    std::cout << "Decision: " << cmdsMap[chosen] << "    Target result: " << *curr_cost << std::endl;

    // Set the exp specific commands
    if (exp_num == 1 && stage % 14 == 0){
        if (target == 0 || target == 1){
            std::string cmd = fmt::format("read_aiger {}; {}; write_aiger {};print_stats;", node.aig_addr.string(), add_aig_opt_cmd, node.aig_addr.string());
            aig_res ress = res_aig(remove_ansi_sequences(executeAbcCommand(cmd)));
            node.aig_node = ress.aig_node;
            node.aig_lev = ress.levs;
            node.seq.push_back(Oper(Cmd::SPECIAL_AIG));

            std::cout<<"    aig  node: "<<node.aig_node<<"    aig level: "<< node.aig_lev << std::endl;
        } else {
            std::string cmd = fmt::format("read_aiger {}; {}; write_aiger {}; if -K 6; print_stats;", node.aig_addr.string(), add_fpga_opt_cmd, node.aig_addr.string());
            map_res ress = res_lut(remove_ansi_sequences(executeAbcCommand(cmd)));
            node.lut_area = ress.lut_nd;
            node.lut_lev = ress.lut_lev;
            node.seq.push_back(Oper(Cmd::SPECIAL_LUT));

            std::cout<<"    lut  area: "<< node.lut_area <<"    lut level: "<< node.lut_lev << std::endl;
        }        
    } else {
        if (target == 0 || target == 1){
            std::cout<<"aig node: "<< node.aig_node<<"    aig level: "<< node.aig_lev << std::endl;
        } else{
            std::cout<<"lut area: "<< node.lut_area<<"    lut level: "<< node.lut_lev << std::endl;
        }
    }
    
    delete(C);
    return chosen; 
}

}