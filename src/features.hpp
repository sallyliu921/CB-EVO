#ifndef CIRCUIT_FEATURES_H
#define CIRCUIT_FEATURES_H

#include <iostream>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <string>
#include <filesystem>
#include <unordered_map>
#include <fstream>
#include <Eigen/Dense>
#include <fmt/core.h>
#include "node.hpp"
#include "utils.hpp"

namespace fs = std::filesystem;

namespace Bandit{

struct TurningPoints {
    int shape_type;
    float range_mean;
    float range_variance;

};

class CircuitFeatures {
public:
    CircuitFeatures() {}

    std::unordered_map<std::string, float> ccirc_stats(const std::string& design_file, const std::string& ccirc_binary, int stage, Cmd arm) {
        std::unordered_map<std::string, float> stats;

        std::string stats_file = fmt::format("./netlist_{}_{}.stats", stage, static_cast<int>(arm));
        std::string ccirc_command = ccirc_binary + " " + design_file + " --partitions 1 --out " + stats_file + " > /dev/null 2>&1";

        // Run ccirc command
        try {
            std::system(ccirc_command.c_str());

            std::string readfile = stats_file;
            std::ifstream file(readfile);
            std::string line;
            while (std::getline(file, line)) {
                if (line.find("Number_of_Nodes") != std::string::npos) {
                    stats["Number_of_Nodes"] = std::stoi(line.substr(line.find_last_of(' ') + 1));          
                }
                if (line.find("Number_of_Edges") != std::string::npos) {
                    stats["Number_of_Edges"] = std::stoi(line.substr(line.find_last_of(' ') + 1));
                }
                if (line.find("Maximum_Delay") != std::string::npos) {
                    stats["Maximum_Delay"] = std::stod(line.substr(line.find_last_of(' ') + 1));
                }
                if (line.find("Number_of_Combinational_Nodes") != std::string::npos) {
                    stats["Number_of_Combinational_Nodes"] = std::stod(line.substr(line.find_last_of(' ') + 1));
                }
                if (line.find("Number_of_high_degree_comb") != std::string::npos) {
                    stats["Number_of_high_degree_comb"] = std::stod(line.substr(line.find_last_of(' ') + 1));
                }
                if (line.find("Number_of_10plus_degree_comb") != std::string::npos) {
                    stats["Number_of_10plus_degree_comb"] = std::stod(line.substr(line.find_last_of(' ') + 1));
                }
                // if (line.find("kin") != std::string::npos) {
                //     stats["kin"] = std::stod(line.substr(line.find_last_of(' ') + 1));
                // }
                if (line.find("Reconvergence:") != std::string::npos) {
                    stats["Reconvergence"] = std::stod(line.substr(line.find_last_of(' ') + 1));
                }
                if (line.find("Reconvergence_max:") != std::string::npos) {
                    stats["Reconvergence_max"] = std::stod(line.substr(line.find_last_of(' ') + 1));    
                }
                if (line.find("Reconvergence_min:") != std::string::npos) {
                    stats["Reconvergence_min"] = std::stod(line.substr(line.find_last_of(' ') + 1));
                }

                if (line.find("Avg_fanout:") != std::string::npos) {
                    extractFloatNumbers(line, "Avg_fanout", stats);
                }
                if (line.find("Avg_fanout_comb:") != std::string::npos) {
                    extractFloatNumbers(line, "Avg_fanout_comb", stats);
                }
                if (line.find("Avg_fanout_pi:") != std::string::npos) {
                    extractFloatNumbers(line, "Avg_fanout_pi", stats);
                }
            
                if (line.find("Node_shape:") != std::string::npos) {
                    process_shape_data(line, "Node_shape", stats);
                }
                if (line.find("Input_shape:") != std::string::npos) {
                    process_shape_data(line, "Input_shape", stats);
                }
                if (line.find("Output_shape:") != std::string::npos) {
                    process_shape_data(line, "Output_shape", stats);
                
                }
                if (line.find("Edge_length_distribution:") != std::string::npos) {
                    process_shape_data(line, "Edge_length_distribution", stats);
                }
                if (line.find("Fanout_distribution:") != std::string::npos) {
                    process_shape_data(line, "Fanout_distribution", stats);
                }

            }
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            return stats;
        }
        
        return stats;
    }

    std::map<std::string, int> yosys_stats(const std::string& design_file, const std::string& yosys_binary) {
        std::map<std::string, int> stats;
        std::string yosys_command = "read_aiger " + design_file + "; stat";
        try {
            std::string command = yosys_binary + " -QT -p '" + yosys_command + "'";
            FILE* pipe = popen(command.c_str(), "r");
            if (!pipe) {
                throw std::runtime_error("Error executing yosys command.");
            }

            char buffer[128];
            std::stringstream output;
            while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
                output << buffer;
            }

            pclose(pipe);

            std::string line;
            while (std::getline(output, line)) {
                if (line.find("Number of wires") != std::string::npos) {
                    stats["number_of_wires"] = std::stoi(line.substr(line.find_last_of(' ') + 1));
                }
                if (line.find("Number of cells") != std::string::npos) {
                    stats["number_of_cells"] = std::stoi(line.substr(line.find_last_of(' ') + 1));
                }
                if (line.find("$_NOT_") != std::string::npos) {
                    stats["nots"] = std::stoi(line.substr(line.find_last_of(' ') + 1));
                }
            }

        } catch (const std::exception& e) {
            std::cout << e.what() << std::endl;
            return std::map<std::string, int>();
        }

        return stats;
    }

    void fix_blif_model_name(const std::string& blif_file) {
        std::ifstream file(blif_file);
        std::vector<std::string> lines;
        std::string line;
        
        // Read all lines
        while (std::getline(file, line)) {
            lines.push_back(line);
        }
        file.close();
        
        // Fix the second line (.model line)
        if (lines.size() >= 2 && lines[1].find(".model ") == 0) {
            std::string model_line = lines[1];
            size_t last_slash = model_line.find_last_of('/');
            if (last_slash != std::string::npos) {
                std::string model_name = model_line.substr(last_slash + 1);
                lines[1] = ".model " + model_name;
            }
        }
        
        // Write back to file
        std::ofstream outfile(blif_file);
        for (const auto& line : lines) {
            outfile << line << '\n';
        }
        outfile.close();
    }

private:

    TurningPoints shape_classifier(const std::vector<float>& datas){
        std::vector<float> gradients;
        for (size_t i = 1;i < datas.size(); i++){
            gradients.push_back(datas[i] - datas[i - 1]);
        }

        int maxima_num = 0;
        std::vector<int> max_locations;
        int count = 0;
        float threshold = 2.0 * calculate_mean(datas);

        for (size_t i = 0; i < gradients.size() - 1; i++){
            count++;
            if ((gradients[i] > threshold / 2.0) && (gradients[count] < threshold / 2.0) 
                    && (gradients[i] != gradients[count])){
                maxima_num++;
                max_locations.push_back(count);
            }
        }

        // Deciding the types
        int shape_type = 0;
        if (maxima_num == 0) {
            shape_type = 0; // Decreasing
        }
        else if (maxima_num <= 5) {
            shape_type = 1; // Cornical
        }
        else {
            shape_type = 2; // Others
        }
        
        // Calculate range data statistics
        float range_mean = std::accumulate(datas.begin(), datas.end(), 0.0) / datas.size();

        float sum_squared_diff = 0.0;
        for (const auto& data : datas) {
            sum_squared_diff += std::pow(data - range_mean, 2);
        }
        float range_variance = sum_squared_diff / datas.size();

        // Add range data statistics to the result
        TurningPoints turning_points = {shape_type, range_mean, range_variance};

        return turning_points;
    }
    
    float calculate_mean(const std::vector<float>& data){
        float sum = 0.0;
        for (const auto& value : data){
            sum += value;
        }
        return sum / data.size();
    }

    void extractFloatNumbers(const std::string& line, const std::string& key, std::unordered_map<std::string, float>& stats) {
        size_t start = line.find("(");
        size_t end = line.find(")");
        std::string insideParentheses = line.substr(start + 1, end - start - 1);

        size_t colonPos = line.find(":");
        std::string outsideParentheses = line.substr(colonPos + 2, start - colonPos - 3);

        float insideValue, outsideValue;
        std::istringstream issInside(insideParentheses);
        std::istringstream issOutside(outsideParentheses);
        issInside >> insideValue;
        issOutside >> outsideValue;

        stats[key + "_1"] = insideValue;
        stats[key + "_2"] = outsideValue;
    }

    void process_shape_data(const std::string& line, const std::string& key_prefix, std::unordered_map<std::string, float>& stats) {
        
        size_t startPos = line.find("(") + 1;
        size_t endPos = line.size() - 1;
        std::string numStr = line.substr(startPos, endPos - startPos);

        std::vector<float> node_shape_data;
        std::istringstream iss(numStr);
        int num;

        while (iss >> num) {
            node_shape_data.push_back(num);
        }


        TurningPoints turning_points = shape_classifier(node_shape_data);

        stats[key_prefix + "_1"] = static_cast<float>(turning_points.shape_type);
        stats[key_prefix + "_2"] = turning_points.range_mean;
        stats[key_prefix + "_3"] = turning_points.range_variance/1000;
    
    }

};

std::vector<float> extractInputData(Node& st_node, int stage, Cmd arm, bool use_ccirc) {
    CircuitFeatures circuit;
    std::vector<float> input_data;

    if (!use_ccirc) {
        // call yosys 
        std::string  yosys_binary = "yosys"; 

        std::map<std::string, int> stats = circuit.yosys_stats(st_node.aig_addr.string(), yosys_binary);
        
        float num_wires = stats["number_of_wires"];      
        float num_cells = stats["number_of_cells"];
        float orig_not = stats["nots"];   

        // prepare for saving into the context vector
        input_data.insert(input_data.end(), {static_cast<float>(num_wires), static_cast<float>(num_cells),
                                            static_cast<float>(orig_not)});
        

    } else {

        // call ccirc
        std::string cmd_blif = fmt::format("read_aiger {}; write_blif netlist_{}.blif", st_node.aig_addr.string(), static_cast<int>(arm));
        std::string output = executeAbcCommand(cmd_blif);
        
        std::string ccircBinary = "./ccirc"; 
        std::string design_file = fmt::format("./netlist_{}.blif", static_cast<int>(arm));
        
        // Fix the .model name in BLIF file
        circuit.fix_blif_model_name(design_file);

        std::unordered_map<std::string, float> stats = circuit.ccirc_stats(design_file, ccircBinary, stage, arm);
        std::vector<std::string> keys = {
            "Number_of_Nodes", "Number_of_Edges", "Maximum_Delay",
            "Number_of_Combinational_Nodes", 
            "Reconvergence", "Reconvergence_max", "Reconvergence_min",
            "Avg_fanout_1", "Avg_fanout_2", "Avg_fanout_comb_1", "Avg_fanout_comb_2", "Avg_fanout_pi_1", "Avg_fanout_pi_2",
            "Node_shape_2","Input_shape_2","Input_shape_3",
            "Output_shape_2", "Output_shape_3", "Edge_length_distribution_2",
            "Edge_length_distribution_3", "Fanout_distribution_2", "Fanout_distribution_3"
        };

        for (const auto& key : keys) {
            input_data.push_back(stats[key]);
        }
    }
    
    return input_data;
}

}
#endif