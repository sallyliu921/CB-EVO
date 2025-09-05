#pragma once
#include <iostream>
#include <sstream>
#include <string>
#include <regex>
#include <utility>    
#include <vector>
#include <climits>
#include <limits>
#include <type_traits>
#include <unordered_map>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <chrono>
#include <cstdlib>
#include "cmd.hpp"

using HashTable = std::unordered_map<int, std::unordered_map<int, int>>;
template<typename T>
class LongTermStatsMap {
private:
    std::unordered_map<int, std::unordered_map<Cmd, T>> data;
    std::unordered_map<int, std::unordered_map<Cmd, float>> scale_data;

public:
    LongTermStatsMap() {
        if constexpr (std::is_same_v<T, int>) {
            // For int types initialized to INT_MAX
            for (auto& pair : data) {
                for (auto& innerPair : pair.second) {
                    innerPair.second = INT_MAX;
                }
            }
        } else if constexpr (std::is_same_v<T, float>) {
            // For float types initialized to 0.0
            for (auto& pair : data) {
                for (auto& innerPair : pair.second) {
                    innerPair.second = 0.0f;
                }
            }
        }
    }

    void update(int id, Cmd command, T value) {
        data[id][command] = value;
    }

    T get(int id, Cmd command) const {
        auto idIt = data.find(id);
        if (idIt != data.end()) {
            auto cmdIt = idIt->second.find(command);
            if (cmdIt != idIt->second.end()) {
                return cmdIt->second;
            }
        }
        throw std::runtime_error("Command not found.");
    }

    bool exists(int id, Cmd command) const {
        auto idIt = data.find(id);
        if (idIt != data.end()) {
            return idIt->second.find(command) != idIt->second.end();
        }
        return false;
    }

    T findMinForId(int id) const {
        auto idIt = data.find(id);
        if (idIt == data.end()) {
            throw std::runtime_error("ID not found.");
        }

        T minVal = std::numeric_limits<T>::max();  
        for (const auto& cmdPair : idIt->second) {
            if (cmdPair.second < minVal) {
                minVal = cmdPair.second;
            }
        }

        return minVal;
    }

    std::pair<float, float> findExtremalForId(int id) const {
        auto idIt = scale_data.find(id);
        if (idIt == scale_data.end()) {
            throw std::runtime_error("ID not found.");
        }

        float maxVal = std::numeric_limits<float>::lowest();  
        float minVal = std::numeric_limits<float>::max();
        for (const auto& cmdPair : idIt->second) {
            if (cmdPair.second > maxVal) {
                maxVal = cmdPair.second;
            }
            if (cmdPair.second < minVal) {
                minVal = cmdPair.second;
            }
        }

        return std::make_pair(maxVal, minVal);
    }

    void scaleData(int id, float scale_ratio) {
        // Get the sum of all the data
        T sum = T();
        bool allSame = true;
        auto idIt = data.find(id);
        std::size_t size;
        if (idIt != data.end()) {
            auto& cmdMap = idIt->second;
            size = cmdMap.size();
            if (!cmdMap.empty()) {
                auto iter = cmdMap.begin();
                T firstValue = iter->second;
                sum += firstValue;  

                for (++iter; iter != cmdMap.end(); ++iter) {
                    sum += iter->second;  
                    if (iter->second != firstValue) {
                        allSame = false;  
                    }
                }
            } else {
                throw std::runtime_error("Data map is empty!");
            }
        } else {
            throw std::runtime_error("Data map is empty!");
        }

        // Scale the data
        if (allSame) {
            // std::cout << "Data fully same!"<<std::endl;
            auto& cmdMap = idIt->second;
            for (const auto& pair : cmdMap) {
                scale_data[id][pair.first] = 1.0f;
            }
        } else {
            float average = static_cast<float>(sum / size);

            // Scale each value based on the average
            auto& cmdMap = idIt->second;
            for (const auto& pair : cmdMap) {
                float scaledValue = average / pair.second;
                scaledValue = pow(scaledValue, scale_ratio);
                scale_data[id][pair.first] = scaledValue;
            }

            auto extremalVal = findExtremalForId(id);
            float range = extremalVal.first - extremalVal.second;
            
            auto& cmdMap_s = scale_data.find(id)->second;
            for (const auto& pair : cmdMap_s) {
                float extendDiffVal = (pair.second - extremalVal.second) / range * 2.0f;
                scale_data[id][pair.first] = extendDiffVal + 0.001f;
            }
        }
    }

    float getScaledData(int id, Cmd command) const {
        auto idIt = scale_data.find(id);
        if (idIt != scale_data.end()) {
            auto cmdIt = idIt->second.find(command);
            if (cmdIt != idIt->second.end()) {
                return cmdIt->second;
            }
        } else{
            throw std::runtime_error("Command not found.");
        }
        return 0;
    }

    void printData() const {
        for (const auto& outer_pair : data) {
            std::cout << "Search depth: " << outer_pair.first << std::endl;
            for (const auto& inner_pair : outer_pair.second) {
                std::cout << "  Arm: " << static_cast<int>(inner_pair.first) << ", Value: " << inner_pair.second << std::endl;
            }
        }
    }
};

struct aig_res{
    int pis;
    int pos;
    int aig_node;
    int levs;
};

struct map_res{
    int aig_node = 0;
    int aig_lev = 0;
    int lut_nd = 0;
    int lut_lev = 0;
};

aig_res res_aig(const std::string& input) {
    // Regular expressions for various parameters
    std::regex re_io(R"(: i/o\s*=\s*(\d+)/\s*(\d+))");
    std::regex re_and(R"(\band\s*=\s*(\d+))");
    std::regex re_lev(R"(\blev\s*=\s*(\d+))");

    std::smatch match;
    int pis = -1, pos = -1, and_nodes = -1, levs = -1;

    // Match and extract the number of primary inputs and outputs
    if (std::regex_search(input, match, re_io) && match.size() > 2) {
        pis = std::stoi(match.str(1));
        pos = std::stoi(match.str(2));
    }

    // Match and extract the number of AND gates
    if (std::regex_search(input, match, re_and) && match.size() > 1) {
        and_nodes = std::stoi(match.str(1));
    }

    // Match and extract the number of levels
    if (std::regex_search(input, match, re_lev) && match.size() > 1) {
        levs = std::stoi(match.str(1));
    }
    return {pis, pos, and_nodes, levs};
}

map_res res_lut(const std::string& input) {
    map_res result;

    std::vector<std::string> lines;
    std::istringstream stream(input);
    std::string line;
    while (std::getline(stream, line)) {
        lines.push_back(line);
    }

    std::regex lev_regex("lev\\s*=\\s*(\\d+)");
    std::regex nd_regex("nd\\s*=\\s*(\\d+)");
    // std::regex and_regex("and\\s*=\\s*(\\d+)");

    std::smatch match;

    for (size_t i = 0; i < lines.size(); ++i) {
        bool aig_flag = false;
        if (std::regex_search(lines[i], match, nd_regex)) {
            result.lut_nd = std::stoi(match[1].str());
        }

        // if (std::regex_search(lines[i], match, and_regex)) {
        //     result.aig_node = std::stoi(match[1].str()); 
        //     aig_flag = 1;
        // }

        if (std::regex_search(lines[i], match, lev_regex)) {
            int lev_value = std::stoi(match[1].str());
            if (aig_flag == 1) {
                result.aig_lev = lev_value;
            } else {
                result.lut_lev = lev_value;
            }
        }
    }

    // std::cout << "input:" <<input<<std::endl;
    // std::cout <<fmt::format("{},{},{},{}", result.aig_node,result.aig_lev,result.lut_nd,result.lut_lev);
    return result;
}

map_res res_lut_new(const std::string& input) {
    map_res result;

    // Regular expressions to extract the needed values
    std::regex and_pattern(R"(\band\s*=\s*(\d+))");
    std::regex aig_lev_pattern(R"(\blev\s*=\s*(\d+))");
    std::regex lut_pattern(R"(\blut\s*=\s*(\d+))");
    std::regex lut_lev_pattern(R"(\blev\s*=\s*(\d+))");

    std::smatch matches;

    // Extract "and" (aig_node)
    if (std::regex_search(input, matches, and_pattern)) {
        result.aig_node = std::stoi(matches[1].str());
    }

    // Extract the first "lev" (aig_lev)
    auto mapping_pos = input.find("Mapping");
    std::string before_mapping = input.substr(0, mapping_pos);
    if (std::regex_search(before_mapping, matches, aig_lev_pattern)) {
        result.aig_lev = std::stoi(matches[1].str());
    }

    // Extract "lut" (lut_node)
    if (std::regex_search(input, matches, lut_pattern)) {
        result.lut_nd = std::stoi(matches[1].str());
    }

    // Extract the second "lev" (lut_lev) after "Mapping" keyword
    std::string after_mapping = input.substr(mapping_pos);
    if (std::regex_search(after_mapping, matches, lut_lev_pattern)) {
        result.lut_lev = std::stoi(matches[1].str());
    }
   
    return result;
}

std::string remove_ansi_sequences(const std::string& input) {
    std::regex ansi_pattern("\033\\[[0-9;]*m");
    return std::regex_replace(input, ansi_pattern, "");
}

std::string extract_case_name(const std::string& input_path) {
    std::regex re(R"([^/]+(?=\.(aig|blif)$))");
    std::smatch match;

    if (std::regex_search(input_path, match, re) && !match.empty()) {
        return match.str(0);  
    }

    return ""; 
}

// init seq file 
struct DataInfo {
    std::string sequence;
    int area;
    int level;
};

struct VectorHash {
    size_t operator()(const std::vector<int>& v) const {
        std::hash<int> hasher;
        size_t seed = 0;
        for (int i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};


// Convert the string sequence to independent Opers
std::string strip(const std::string &s){
    std::string res;
    for(int i = 0,valid = 0; i < s.size(); ++i){
        if(!isspace(s[i]))
            valid=1;
        if(valid)
            res.push_back(s[i]);
    }
    while(!res.empty() && isspace(res.back()))
        res.pop_back();
    return res;
}

std::vector<std::string> split(const std::string &s,char delimiter){
    std::vector<std::string> tokens;
    std::string token;
    for(char c:s){
        if(c == delimiter){
            tokens.emplace_back(strip(token));
            token.clear();
        }
        else token.push_back(c);
    }
    if(!token.empty())
        tokens.emplace_back(strip(token));
    return tokens;
}

int cnt_seq_len(const std::string& seq){
    int count = 0;
    
    for (char ch : seq) {
        if (ch == ';') {
            count++;
        }
    }
    return count;
    
}

std::string executeAbcCommand(const std::string& cmd) {
    // Create the full command string
    std::string fullCommand = "./abc -c \"" + cmd + "\"";
    std::string result;

    try {
        // Use unique_ptr to manage the FILE* from popen and ensure pclose is called
        std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(fullCommand.c_str(), "r"), pclose);
        if (!pipe) {
            throw std::runtime_error("popen() failed!");
        }

        // Read output to a string
        char buffer[256];
        while (fgets(buffer, sizeof(buffer), pipe.get()) != nullptr) {
            result += buffer;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return ""; // Return empty string on error
    }

    return result;
}

/**
 * @brief the tic_toc class mainly used to calcute the procedure's execute time.
 */
class tic_toc{
public:
    tic_toc()
    {
    tic();
    }

    inline void tic()
    {
    _start = std::chrono::system_clock::now();
    }

    inline double toc()
    {
    _end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = _end - _start;
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(elapsed_seconds);
    return double( duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
    }
private:
    std::chrono::time_point<std::chrono::system_clock> _start , _end;
};  

std::string removeSpaces(const std::string& input) {
    std::string output;
    for (char c : input) {
        if (c != ' ') {
            output.push_back(c);
        }
    }
    return output;
}

std::pair<float, float> findExtremalForId(const std::unordered_map<Cmd, float>& cost) {
    float maxVal = std::numeric_limits<float>::lowest();  
    float minVal = std::numeric_limits<float>::max();

    // Iterate over the map to find max and min
    for (const auto& cmdPair : cost) {
        if (cmdPair.second > maxVal) {
            maxVal = cmdPair.second;
        }
        if (cmdPair.second < minVal) {
            minVal = cmdPair.second;
        }
    }

    return std::make_pair(maxVal, minVal);
}

void printHashTable(const HashTable& valuetable, int curr_res) {
    std::cout << "Printing valuetable:" << std::endl;
    for (const auto& kv : valuetable) {
        int stage = kv.first;
        const std::unordered_map<int, int>& innerMap = kv.second;
        std::cout << "stage: " << stage << std::endl;
        for (const auto& innerKv : innerMap) {
            int stage = innerKv.first;
            double value = innerKv.second;
            std::cout << "Search_depth: " << stage << ", Value: " << value << std::endl;
        }
    }
    std::cout << "Printing curr_res:" << curr_res << std::endl;
}

uint32_t rand_int(uint32_t l,uint32_t r){
    return rand()%(r-l+1)+l;
}
float rand_float(float l, float r) {
    uint32_t precision = 1000000;
    float random_fraction = static_cast<float>(rand_int(0, precision)) / precision;
    return l + (r - l) * random_fraction;
}
