#include "Parser.h"
#include "FunctionRealization.h"

int main(int argc, char * argv[])
{
    Parser parser = Parser(argc, argv);
    FunctionParameters params = parser.Parse();
    FunctionRealization solve = FunctionRealization(params);
    time_t start, end;
    time(&start);
    solve.ComputeFunction();
    time(&end);
    double time_taken = double(end - start); 
    cout << "Time taken by program is : " << fixed 
        << time_taken << " sec " << endl;
    solve.SaveFunctionsValuesToTxtFiles();
    return 0;
}