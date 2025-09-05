#pragma once
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <stdexcept>

enum class Cmd {
    RW, RWZ, B, RF, RFZ, RESUB, RESUBZ, NONE, SPECIAL_AIG, SPECIAL_LUT, CHOICE_S, CHOICE_C, PUT
};

inline std::ostream& operator<<(std::ostream& os, const Cmd& cmd) {
    static const std::map<Cmd, std::string> cmdNames = {
        {Cmd::RW, "rewrite"},
        {Cmd::RWZ, "rewrite -z"},
        {Cmd::B, "balance"},
        {Cmd::RF, "refactor"},
        {Cmd::RFZ, "refactor -z"},
        {Cmd::RESUB, "resub"},
        {Cmd::RESUBZ, "resub -z"},
        {Cmd::NONE, "NONE"},
        {Cmd::SPECIAL_AIG, "SPECIAL_AIG"},
        {Cmd::SPECIAL_LUT, "SPECIAL_LUT"},
        {Cmd::CHOICE_S, "&choice_store"},
        {Cmd::CHOICE_C, "&choice_compute"},
        {Cmd::PUT, "&put"}
    };
    
    auto it = cmdNames.find(cmd);
    if (it != cmdNames.end()) {
        os << it->second;
    } else {
        os << "UNKNOWN";
    }
    return os;
}

class CmdManager {
private:
    std::map<Cmd, std::string> typeToStr;
    std::vector<Cmd> arms;

    static std::map<Cmd, std::string> getDefaultCommands() {
        return {
            {Cmd::RW, "rewrite"},
            {Cmd::RWZ, "rewrite -z"},
            {Cmd::B, "balance"},
            {Cmd::RF, "refactor"},
            {Cmd::RFZ, "refactor -z"},
            {Cmd::RESUB, "resub"},
            {Cmd::RESUBZ, "resub -z"}
        };
    }

public:
    CmdManager() {
        typeToStr = getDefaultCommands();
        updateArms();  
    }

    CmdManager(const CmdManager& other) : typeToStr(other.typeToStr), arms(other.arms) {}
    
    CmdManager& operator=(const CmdManager& other) {
        if (this != &other) {
            typeToStr = other.typeToStr;
            arms = other.arms;
        }
        return *this;
    }

    void updateArms() {
        arms.clear(); 
        for (const auto& cmdPair : typeToStr) {
            arms.push_back(cmdPair.first);
        }
    }

    const std::vector<Cmd>& getArms() const {
        return arms;
    }

    std::map<Cmd, std::string> getCurrentTypeToStr() const {
        return typeToStr;
    }

    int getCurrentMapSize() const {
        return typeToStr.size();
    }

    std::string getFullCommand(Cmd cmd) {
        auto currentCmds = getCurrentTypeToStr();
        auto it = currentCmds.find(cmd);
        if (it != currentCmds.end()) {
            return it->second;
        } else if (cmd == Cmd::SPECIAL_AIG || cmd == Cmd::SPECIAL_LUT){
            return "";
        } else if (cmd == Cmd::CHOICE_S){
            return "&get; &choice_store";
        } else if (cmd == Cmd::CHOICE_C){
            return "&get; &choice_compute";
        } else if (cmd == Cmd::PUT){
            return "&put";
        } else if (cmd == Cmd::NONE){
            return "";
        } else{
            throw std::runtime_error("Command has been removed or does not exist.");
        }
    }

    int getCmdIndex(Cmd cmd) const {
        auto currentCmds = getCurrentTypeToStr(); 
        int index = 0; 
        for (const auto& pair : currentCmds) 
        {
            if (pair.first == cmd)
                return index;
            index++; 
        }
        throw std::runtime_error("Command does not exist or has been removed.");
    }

    void removeCmd(Cmd cmd) {
        // First collect all commands and their strings in order
        std::vector<std::pair<Cmd, std::string>> cmds;
        for (const auto& pair : typeToStr) {
            if (pair.first != cmd) {
                cmds.push_back(pair);
            }
        }

        // Clear and rebuild typeToStr with continuous indices
        typeToStr.clear();
        for (size_t i = 0; i < cmds.size(); i++) {
            typeToStr[static_cast<Cmd>(i)] = std::move(cmds[i].second);
        }
        
        // Update arms vector to match typeToStr
        updateArms();
    }
};

