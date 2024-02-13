
#include "../include/Solver.hpp"

int main(int arcg, char *argv[])
{
    Scalar L0 = 10;
    Scalar N = 5;
    Scalar D_outer = 0.1;
    Scalar D_inner = 0.05;

    Config config;
    config.n = 2;
    config.dt = 0.01;
    config.t = 0;

    Geometry geometry{L0, N, D_outer, D_inner};

    cout << "hello\n";
    solve(config, geometry);
}