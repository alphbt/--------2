#include "Parser.h"
#include "FunctionRealization.h"
#include <iomanip>

int main(int argc, char * argv[])
{
    Parser parser = Parser(argc, argv);
    FunctionParameters params = parser.Parse();
    FunctionRealization solve = FunctionRealization(params);
    clock_t start, end;
    start = clock();
    solve.ComputeFunction();
    end = clock();
    double time_taken = double(end - start) / CLOCKS_PER_SEC; 
    cout << "Time taken by program is : " << fixed 
        << time_taken << setprecision(5) << " sec " << endl;
    solve.SaveFunctionsValuesToTxtFiles();
    return 0;
}