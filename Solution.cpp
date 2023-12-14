#include "Parser.h"
#include "FunctionRealization.h"
#include <iomanip>

int main(int argc, char * argv[])
{
    Parser parser = Parser(argc, argv);
    FunctionParameters params = parser.Parse();
    FunctionRealization solve = FunctionRealization(params);
    double start, end;
    start = omp_get_wtime();
    solve.ComputeFunction();
    end = omp_get_wtime();
    setprecision(9);
    double time_taken = end - start; 
    cout << "Time taken by program is : " << time_taken << " sec " << endl; 
    double err = solve.GetMaxError();
    cout << "Max error is : " << err << endl;
    solve.SaveFunctionsValuesToTxtFiles();
    return 0;
}