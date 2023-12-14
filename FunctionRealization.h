#include "omp.h"
#include "FunctionParameters.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;


class FunctionRealization{
    private:
        FunctionParameters params = FunctionParameters();
        vector<vector<vector<vector<long double>>>> u;
        vector<double> error;
        double a2 = 1 / (4 * M_PI * M_PI);
        void SetBordersValues();
        void SetInteriorValues();
        void SetZeroTimeValues();
        void SetOneTimeValues();
        double GetAnalyticalSolve(double x, double y, double z, double t);
        double GetLaplasian(int i, int j, int k, int t);
        void SaveErrorFile();
        void SaveProjectionsFiles();
        void ComputeError();

    public:
        FunctionRealization(const FunctionParameters &params);
        void ComputeFunction();
        double GetMaxError();
        void SaveFunctionsValuesToTxtFiles();
};