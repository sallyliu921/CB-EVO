#include <filesystem>
#include <iostream>
#include <climits>
#include <time.h>
#include <algorithm>
#include <CLI/CLI.hpp>
#include <fmt/core.h>
#include "cmd.hpp"
#include "utils.hpp"
#include "bandit.hpp"
#include "evolution.hpp"

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    CLI::App app{"CB-EVO generator"};

    // Global variables defined in `node.hpp`
    app.add_option("-i", input_file, "Input file path")
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("-t", target, "Target: 0=AIG nodes, 1=AIG nodes*level, 2=LUTs, 3=LUTs*level")
        ->required()
        ->check(CLI::Range(0, 3));
    app.add_option("-s", seq_all, "Total sequence length")
        ->check(CLI::Range(1, 400));
    app.add_option("-e", exp_num, "Experiment number")
        ->check(CLI::Range(0, 1));
    app.add_flag("-f", use_static_extraction, "Enable ccirc for additional static feature extraction in bandit model (default: false)");

    CLI11_PARSE(app, argc, argv);
    std::string case_name = extract_case_name(input_file);  
    if (use_static_extraction && (std::find(vtr_cases.begin(), vtr_cases.end(), case_name) != vtr_cases.end())) {
        exp_num = 1;
    }
    std::cout << "Current case name: " << case_name << std::endl;

    if (seq_all != 0){
        std::cout << "Using fixed sequence length: " << seq_all << std::endl;
    } else {
        std::cout << "Using dynamic sequence length: Bandit disabled, Evo will extend sequence based on case size" << std::endl;
    }

    // tmp folder for saving aig outputs
    fs::remove_all(tmp_dir_path); 
    fs::path path_in = input_file;  
    std::string path_out = "./" + case_name + ".seq"; 
    Node ori_node = get_ori_node(path_in);
    CmdManager manager;
    tic_toc t_clock;

    //////////////// Contextual Bandit Tuning ////////////////////
    Node stage_node = ori_node;

    if (seq_all != 0) {
        // Initialize parameters for Bandit algorithm
        if (exp_num == 1) {
            Bandit::seq_bandit_length = 14;
        } else {
            Bandit::seq_bandit_length = static_cast<int>(seq_all / 4);
        }
        assert(Bandit::seq_bandit_length != 0 && "Error: seq_bandit_length cannot be 0. Please check the value of seq_all.");
                
        // Return-back parameters
        std::vector<Bandit::StageInfo> stage_history(Bandit::seq_bandit_length);  // Attention! Start from 0!
        double rt_threshold = 2.0;
        HashTable valuetable;
        
        // Run contextual bandit tuning
        std::cout << "Starting contextual bandit tuning!"<< std::endl;
        std::cout << "******************************************************"<<std::endl;
        
        for (int stage = 1; stage <= Bandit::seq_bandit_length; stage++){
            std::cout << fmt::format("Stage {}{} is working...", stage, (stage_history[stage-1].needs_retry ? "(rt)" : "")) << std::endl;

            CmdManager tmp_manager = manager;
            if (stage_history[stage-1].needs_retry){
                std::cout << "Removing command: " << stage_history[stage-1].cmd << std::endl;
                tmp_manager.removeCmd(stage_history[stage-1].cmd);
            }    

            int curr_cost = 0;
            Cmd res_cmd = Bandit::Agent(stage_node, tmp_manager, stage, valuetable, &curr_cost);
            stage_history[stage-1].cmd = res_cmd;
            stage_history[stage-1].node = stage_node;
            write_ans(stage_node, path_out, manager);

            // Return back strategy: check if we need to return to a previous stage
            if (!stage_history[stage-1].needs_retry && stage > 1){
                int min_value = INT_MAX; 
                int return_back_stage = -1; // The stage index to return back to

                // Search for minimum value
                for (const auto& [s, innerMap] : valuetable) {
                    // Skip stages that have already been revisited
                    if (s > 0 && stage_history[s-1].needs_retry) {
                        continue;
                    }
                    
                    if (auto it = innerMap.find(stage-1); it != innerMap.end()) {
                        if (it->second < min_value) {
                            min_value = it->second;
                            return_back_stage = s;
                        }
                    }
                }

                std::cout << "\n=== Return Back Mechanism Debug Info ===" << std::endl;
                double diff_ratio = return_back_stage != -1 && min_value != 0 ? (static_cast<double>(curr_cost - min_value) / min_value) : 0.0;
                std::cout << fmt::format("Current Stage: {}, Current Cost: {}, Min Value Stage: {}, Min Value: {}, Diff Ratio: {:.2f}%", 
                                         stage, curr_cost, return_back_stage != -1 ? return_back_stage : -1, min_value, diff_ratio * 100.0) << std::endl;
                assert(return_back_stage < stage && "Error: Cannot return to a stage that hasn't been processed yet");

                // Search must be restarted from the stage following min_s (min_s >= 1).
                if (return_back_stage >= 1 && !stage_history[return_back_stage].needs_retry && (diff_ratio * 100.0 > rt_threshold)){ 
                    // Go back to stage return_back_stage + 1 (index return_back_stage in stage_history)
                    std::cout << "!!! Return back to stage:" << return_back_stage + 1 << std::endl;
                    
                    stage_node = stage_history[return_back_stage].node; 
                    fs::path curr_aig_addr = stage_node.aig_addr;
                    fs::path children_dir = curr_aig_addr.parent_path() / "children";
                    if (fs::exists(children_dir)) { 
                        fs::remove_all(children_dir);   
                    } 
                    
                    stage_history[return_back_stage].needs_retry = true;
                    stage = return_back_stage;    // Increment by 1 at the next loop start.
                } 
            }

            std::cout << "******************************************************"<<std::endl;
        }
        std::cout << "Finish contextual bandit tuning and start evo tuning!" << std::endl;
    } else {
        std::cout << "Skipping contextual bandit tuning (seq_all = 0)" << std::endl;
        std::cout << "******************************************************"<<std::endl;
    }

    //////////////// Evolution algorithm ////////////////////
    // Start searching the optimized sequence
    std::cout << "Starting evolution tuning!"<< std::endl;
    std::cout << "******************************************************"<<std::endl;
    uint32_t limit_population = 6;
    uint32_t select_rate = 7;

    Node evo_node = Evolution::evo_search(manager, stage_node, path_out, limit_population, select_rate);  
    std::cout << "Finish the evolution algorithm!" << std::endl;
    std::cout << "******************************************************" << std::endl;

    //////////////// Final result /////////////////////////
    std::cout << fmt::format("Best Flow decided by CB-EVO: {}", evo_node.get_seq_str(manager)) << std::endl;
    if (target == 0 || target == 1){
        std::cout<<"aig node: "<< evo_node.aig_node<<"    aig level: "<< evo_node.aig_lev << std::endl;
    } else{
        std::cout<<"lut area: "<< evo_node.lut_area<<"    lut level: "<< evo_node.lut_lev << std::endl;
    }
    std::cout<<"Total decision and optimization time : "<<t_clock.toc()<<"s"<<std::endl;      

    fs::copy(evo_node.aig_addr,  fs::current_path() / fmt::format("{}_cb_evo.aiger", case_name), fs::copy_options::overwrite_existing);
    fs::remove_all(tmp_dir_path); 

    return 0;
}