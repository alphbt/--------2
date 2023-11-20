#include "Parser.h"
#include <cmath>

Parser::Parser(int argc, char **argv){
    for (int i=1; i < argc; ++i)
        this->tokens.push_back(std::string(argv[i]));
}

const std::string&  Parser::GetCmdOption(const std::string &option) const{
    std::vector<std::string>::const_iterator itr;
    itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
    if (itr != this->tokens.end() && ++itr != this->tokens.end()){
        return *itr;
    }
    static const std::string empty_string("");
    return empty_string;
}

FunctionParameters Parser::Parse(){
    FunctionParameters params = FunctionParameters();
    params.N = int_parse(GetCmdOption("N"), 50);
    params.T = int_parse(GetCmdOption("T"), 1);
    params.K = int_parse(GetCmdOption("K"), 100);
    params.Lx = double_parse(GetCmdOption("Lx"), 1);
    params.Ly = double_parse(GetCmdOption("Ly"), 1);
    params.Lz = double_parse(GetCmdOption("Lz"), 1);
    params.ComputeDependParameters();
    return params;
} 

int Parser::int_parse(const std::string &option, int default_value){
    if (option == "") return default_value;
    try{
        return std::stoi(option);
    }
    catch(...){
        return default_value;
    }
}

long double Parser::double_parse(const std::string &option, long double default_value){
    if (option == "") return default_value;
    if (option == "pi") return M_PI;
    try{
        return std::stold(option);
    }
    catch(...){
        return default_value;
    }
}