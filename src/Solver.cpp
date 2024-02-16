#include "../include/Solver.hpp"
#include "../include/SolverRuntime.hpp"

void solve(Config &config, const Geometry &geometry)
{
    create_output_dir(config);
    BeamSystem beam_system{config, geometry};
    calc_dt(config, geometry);
    Index n_steps = config.get_n_steps();

    calc_static_loads(config, geometry, beam_system.R_static);

    cout << "\n-----------------Starting simulation----------------\n";
    cout << "Running with dt = " << config.dt << ", for " << n_steps << " timesteps\n"
         << "for a total time of " << config.dt * n_steps << " seconds\n"
         << "--------------------------------------------------------\n";

    Timer timer;
    timer.start_counter();

    for (Index n = 0; n < n_steps; n++)
    {
        config.t = n * config.dt;
        config.n = n;

        // if (n % 1000 == 0)
        //{
        cout << "n = " << n << ", t = " << config.t << "\n";
        //}

        // save output
        save_csv(config, geometry, beam_system);

        // calculate internal loads
        assemble(config, geometry, beam_system);

        // calculate external loads

        // step central differences (Includes updating displacements and nodal triads/Quaternions)

        set_simple_bc(config.bc_case, geometry, beam_system);
    }
    timer.stop_counter();
    timer.print_elapsed_time();
}

void calc_dt(Config &config, const Geometry &geometry)
{
    Scalar dx_min = std::numeric_limits<Scalar>::max();
    for (Index ie = 0; ie < geometry.get_Ne(); ie++)
    {
        dx_min = min(dx_min, geometry.dx_e(ie));
    }
    assert(dx_min > 0 && dx_min < 100);     // just a check
    Scalar c = sqrt(config.E / config.rho); /*Speed of sound*/
    Scalar CFL = config.CFL;
    assert(CFL < 1 && CFL > 0);
    config.dt = CFL * dx_min / c;
}