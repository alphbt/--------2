#include <bits/stdc++.h>
#include <fstream>

using namespace std;

long double Lx, Ly, Lz;

long double u_analytic(const long double &x, const long double &y, const long double &z, const long double &t){
    return sin(2 * M_PI / Lx * x) * sin(M_PI / Ly * y) * sin(M_PI / Lz * z) * cos(t * sqrt(4/(Lx*Lx) + 1/(Ly*Ly) + 1/(Lz*Lz)) / 2.0 );
}

long double delta_h(const long double &h, const vector<vector<vector<vector<long double>>>> &u, const int& i, const int& j, const int& k, const int& t){
    long double ux, uy, uz;
    ux = u[i-1][j][k][t] - 2*u[i][j][k][t] + u[i+1][j][k][t];
    uy = u[i][j-1][k][t] - 2*u[i][j][k][t] + u[i][j+1][k][t];
    uz = u[i][j][k-1][t] - 2*u[i][j][k][t] + u[i][j][k +1][t];
    return (ux+uy+uz) / (h*h);
} 


int main(int argc, char * argv[])
{
    int N = 50;
    int T = 1;
    int K = 150;
    Lx = 1;
    Ly = 1;
    Lz = 1;
    

    long double hx = Lx / (N*1.0);
    long double hy = Ly / (N*1.0);
    long double hz = Lz / (N*1.0);
    long double h  = min(hx, min(hy, hz));

    long double t0 = (double)T / K;
    

    vector<vector<vector<vector<long double>>>> u(N, vector<vector<vector<long double>>>(N, vector<vector<long double>>(N, vector<long double>(K))));
    
    for(int i = 0; i < N; i++)
    {
        for(int k = 0; k < N; k++)
        {
            for(int t = 0; t < K; t++)
            {
                u[i][0][k][t] = 0;
                u[i][N - 1][k][t] = 0;
                u[i][k][0][t] = 0;
                u[i][k][N - 1][t] = 0;
                u[0][i][k][t] = u_analytic(0, i*hy, k*hz, t0*t);
                u[N-1][i][k][t] = u_analytic(Lx, i*hy, k*hz, t0*t);
            }
        }
    }

    for(int i = 1; i < N - 1; i++) // t = 0 u^0
    {
        for(int j = 1; j < N - 1; j++)
        {
            for(int k = 1; k < N - 1; k++)
            {
                u[i][j][k][0] = u_analytic(i*hx, j*hy, k*hz, 0);
            }
        }
    }
    
    long double a2 = 1 / (4 * M_PI * M_PI);

    for(int i = 1; i < N - 1; i++) //t=1 u^1
    {
        for(int j = 1; j < N - 1; j++)
        {
            for(int k = 1; k < N - 1; k++)
            {
                u[i][j][k][1] = u[i][j][k][0] + t0*t0*a2*delta_h(h, u, i, j, k, 0) / 2;
            }
        }
    }

    for(int i = 1; i < N - 1; i++) //t=1 u^1
    {
        for(int j = 1; j < N - 1; j++)
        {
            for(int k = 1; k < N - 1; k++)
            {
                for(int t = 2; t < K; t++)
                {
                    u[i][j][k][t] = 2* u[i][j][k][t-1] - u[i][j][k][t-2] + t0*t0*a2*delta_h(h, u, i, j, k, t-1);
                }
            }
        }
    }

    ofstream errorFile("error.txt");

    for(int t = 0; t < K; t++){
        long double max_error = 0;
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                for(int k = 0; k < N; k++){
                    long double err = abs(u[i][j][k][t] - u_analytic(i*hx, j*hy, k*hx, t0*t));
                    if(err > max_error){
                        max_error = err;
                    } 
                }
            }
        }
        errorFile << max_error << endl;
    }

    ofstream xyFile("xy_data.txt");
    ofstream yzFile("yz_data.txt");
    ofstream xzFile("xz_data.txt");

    vector<long double> xy(N * N);
    int cnt = 0;
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++)
        {
           xyFile << u[i][j][2][2] << endl;
           yzFile << u[2][i][j][2] << endl;
           xzFile << u[i][2][j][2] << endl;
        }
    }
    
    
}