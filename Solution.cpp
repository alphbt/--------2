#include "Parser.h"
#include "FunctionRealization.h"

int main(int argc, char * argv[])
{
    Parser parser = Parser(argc, argv);
    FunctionParameters params = parser.Parse();
    FunctionRealization solve = FunctionRealization(params);
    solve.ComputeFunction();
    solve.SaveFunctionsValuesToTxtFiles();
    return 0;
}