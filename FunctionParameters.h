#pragma once
class FunctionParameters{
    public:
        int N = 50;
        int T = 1;
        int K = 100;
        long double Lx = 1;
        long double Ly = 1;
        long double Lz = 1;
        long double hx;
        long double hy;
        long double hz;
        long double t0;
    
    public:
        FunctionParameters();
        void ComputeDependParameters();
        
};