#include "FunctionParameters.h"

void FunctionParameters::ComputeDependParameters(){
    hx = Lx / N;
    hy = Ly / N;
    hz = Lz / N;
    t0 = T / K;
}

FunctionParameters::FunctionParameters(){
    ComputeDependParameters();
}