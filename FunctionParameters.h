#pragma once
#include <cmath>
class FunctionParameters{
    public:
        int N = 128;
        double T = 0.000001;
        int K = 20;
        double Lx = 1;
        double Ly = 1;
        double Lz = 1;
        double hx;
        double hy;
        double hz;
        double t0;
    
    public:
        FunctionParameters();
        void ComputeDependParameters();
        
};