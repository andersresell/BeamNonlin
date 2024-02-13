#include "../include/Solver.hpp"
#include "../include/SolverRuntime.hpp"

void solve(Config &config, const Geometry &geometry)
{

    BeamSol beam_sol(config.N);

    Timer timer;
    timer.start_counter();

    for (Index n = 0; n < config.n_timesteps; n++)
    {
        config.t = n * config.dt;
        config.n = n;

        // save output

        // calculate internal loads

        // calculate external loads

        // step central differences (Includes updating displacements and nodal triads/quaternions)
    }
    timer.stop_counter();
    timer.print_elapsed_time();
}