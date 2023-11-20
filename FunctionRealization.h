#include "omp.h"
#include "FunctionParameters.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;


class FunctionRealization{
    private:
        FunctionParameters params;
        vector<vector<vector<vector<long double>>>> u;
        long double a2 = 1 / (4 * M_PI * M_PI);
        void SetBordersValues();
        void SetInteriorValues();
        void SetZeroTimeValues();
        void SetOneTimeValues();
        long double GetAnalyticalSolve(long double x, long double y, long double z, long double t);
        long double GetLaplasian(int i, int j, int k, int t);
        void SaveErrorFile();
        void SaveProjectionsFiles();

    public:
        FunctionRealization(const FunctionParameters &params);
        void ComputeFunction();
        void SaveFunctionsValuesToTxtFiles();
};