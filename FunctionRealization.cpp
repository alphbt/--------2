#include "FunctionRealization.h"

FunctionRealization::FunctionRealization(const FunctionParameters &params){ //инициализация параметров и векторов
    this->params.N = params.N;
    this->params.T = params.T;
    this->params.K = params.K;
    this->params.Lx = params.Lx;
    this->params.Ly = params.Ly;
    this->params.Lz = params.Lz;
    this->params.hx = params.hx;
    this->params.hy = params.hy;
    this->params.hz = params.hz;
    this->params.t0 = params.t0;
    u.resize(params.N);
    for(int i = 0; i < params.N; i++){
        u[i].resize(params.N);
        for(int j = 0; j < params.N; j++){
            u[i][j].resize(params.N);
            for(int k = 0; k < params.N; k++){
                u[i][j][k].resize(params.K);
            }
        }
    }
    error.resize(params.K);
}

double FunctionRealization::GetAnalyticalSolve(double x, double y, double z, double t){ // аналитическое решение
    return sin(2 * M_PI / params.Lx * x) * sin(M_PI / params.Ly * y) * sin(M_PI / params.Lz * z) * 
            cos(t * sqrt(4/(params.Lx*params.Lx) + 1/(params.Ly*params.Ly) + 1/(params.Lz*params.Lz)) / 2.0 );
}

void FunctionRealization::SetBordersValues(){ // граничные условия
    #pragma omp parallel for
    for(int i = 0; i < params.N; i++){
        #pragma omp parallel for
        for(int k = 0; k < params.N; k++){
            #pragma omp parallel for
            for(int t = 0; t < params.K; t++){
                u[i][0][k][t] = 0;
                u[i][params.N - 1][k][t] = 0;
                u[i][k][0][t] = 0;
                u[i][k][params.N - 1][t] = 0;
                u[0][i][k][t] = GetAnalyticalSolve(0, i*params.hy, k*params.hz, params.t0*t);
                u[params.N-1][i][k][t] = GetAnalyticalSolve(params.Lx, i*params.hy, k*params.hz, params.t0*t);
            }
        }
    }
}

double FunctionRealization::GetLaplasian(int i, int j, int k, int t){ //лапасиан
    long double ux, uy, uz;
    ux = u[i-1][j][k][t] - 2*u[i][j][k][t] + u[i+1][j][k][t];
    uy = u[i][j-1][k][t] - 2*u[i][j][k][t] + u[i][j+1][k][t];
    uz = u[i][j][k-1][t] - 2*u[i][j][k][t] + u[i][j][k+1][t];
    return ux / params.hx + uy / params.hy + uz / params.hz;
}

void FunctionRealization::SetZeroTimeValues(){ // вычисление u^0
    #pragma omp parallel for
    for(int i = 1; i < params.N - 1; i++){
        #pragma omp parallel for
        for(int j = 1; j < params.N - 1; j++){
            #pragma omp parallel for
            for(int k = 1; k < params.N - 1; k++){
                u[i][j][k][0] = GetAnalyticalSolve(i*params.hx, j*params.hy, k*params.hz, 0);
            }
        }
    }
}

void FunctionRealization::SetOneTimeValues(){ // вычисление u^1
    #pragma omp parallel for collapse(3)
    for(int i = 1; i < params.N - 1; i++){
        for(int j = 1; j < params.N - 1; j++){
            for(int k = 1; k < params.N - 1; k++){
                u[i][j][k][1] = u[i][j][k][0] + params.t0*params.t0*a2*GetLaplasian(i, j, k, 0) / 2;
            }
        }
    }
}

void FunctionRealization::SetInteriorValues(){  // вычисление функции на внутренней сетке
    #pragma omp parallel for collapse(4)
    for(int i = 1; i < params.N - 1; i++){
        for(int j = 1; j < params.N - 1; j++){
            for(int k = 1; k < params.N - 1; k++){
                for(int t = 2; t < params.K; t++){
                    u[i][j][k][t] = 2* u[i][j][k][t-1] - u[i][j][k][t-2] + 
                        params.t0*params.t0*a2*GetLaplasian(i, j, k, t-1);
                }
            }
        }
    }
}

void FunctionRealization::ComputeFunction(){ // основные вычисления
    SetZeroTimeValues();
    SetOneTimeValues();
    SetInteriorValues();
    SetBordersValues();
    ComputeError();
}

void FunctionRealization::ComputeError(){ // подсчет ошибки
    #pragma omp parallel for
    for(int t = 0; t < params.K; t++){
        double max_error = 0;
        #pragma omp parallel for reduction(max:max_error)
        for(int i = 1; i < params.N - 1; i++){
            #pragma omp parallel for reduction(max:max_error)
            for(int j = 1; j < params.N - 1; j++){
                #pragma omp parallel for reduction(max:max_error)
                for(int k = 1; k < params.N - 1; k++){
                    double err = abs(u[i][j][k][t] - GetAnalyticalSolve(i*params.hx, j*params.hy, k*params.hx, params.t0*t));
                    max_error = max(err, max_error);
                }
            }
        }
        error[t] = max_error;
    }
}

double FunctionRealization::GetMaxError(){
    double max_err = error[0];
    for(int i = 1; i < params.K; i++){
        max_err = max(max_err, error[i]);
    }
    return max_err;
}

void FunctionRealization::SaveErrorFile(){
    ofstream errorFile("error.txt");
    
    for(int i = 0; i < params.K; i++){
        errorFile << error[i] << endl;
    }
}

void FunctionRealization::SaveProjectionsFiles(){
    ofstream xyFile("xy_data.txt");
    ofstream yzFile("yz_data.txt");
    ofstream xzFile("xz_data.txt");

    for(int i = 0; i < params.N; i++){
        for(int j = 0; j < params.N; j++)
        {
           xyFile << u[i][j][2][2] << endl;
           yzFile << u[2][i][j][2] << endl;
           xzFile << u[i][2][j][2] << endl;
        }
    }
}

void FunctionRealization::SaveFunctionsValuesToTxtFiles(){
    SaveErrorFile();
    SaveProjectionsFiles();
}