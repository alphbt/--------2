#pragma once
#include "FunctionParameters.h"
#include <algorithm>
#include <string>
#include <vector>

class Parser{
    public:
        Parser(int argc, char **argv);
        const std::string& GetCmdOption(const std::string &option)const;
        FunctionParameters Parse();
        
    private:
        std::vector<std::string> tokens;

        int int_parse(const std::string &option, int default_value);
        long double double_parse(const std::string &option, long double default_value);
};