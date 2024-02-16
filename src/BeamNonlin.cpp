
#include "../include/Solver.hpp"

int main(int arcg, char *argv[])
{

    Scalar L0 = 10;
    Scalar N = 5;
    Scalar D_outer = 0.1;
    Scalar D_inner = 0.05;

    Config config{};
    config.n_max = 20;
    config.t_max = 1;
    config.dt = 0.01;
    config.nu = 0.3;
    config.E = 200 * pow(10, 9);
    config.rho = 7000;
    config.CFL = 0.9;

    Geometry geometry{L0, N, D_outer, D_inner};

    cout << "hello\n";
    solve(config, geometry);
}
