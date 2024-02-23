#include "../include/Solver.hpp"

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

    for (Index n = 0; n <= n_steps; n++)
    {
        config.t = n * config.dt;
        config.n = n;

        // remove later
        n_glob = n;

        if (n % 1000 == 0)
        {
            cout << "\n---------------------------------\n"
                 << "n = " << n << ", t = " << config.t
                 << "\n---------------------------------\n";
        }

        // save output
        save_csv(config, geometry, beam_system);

        // bool set_disp = false;

        // if (set_disp)
        // {
        //     Scalar u_new = -1.0 * (n + 1) * 0.3;
        //     assert(geometry.get_N() == 2);
        //     // beam_system.u[1].trans = {0, u_new, 0};

        //     Scalar theta = 1 * M_PI / 180; // 1 deg
        //     Scalar c = cos(theta);
        //     Scalar s = sin(theta);
        //     Mat3 T;
        //     T << c, 0, -s,
        //         0, 1, 0,
        //         s, 0, c;
        //     Mat3 U = Mat3::Identity();

        //     beam_system.u[0].rot.from_matrix(T);
        //     beam_system.u[1].rot.from_matrix(U);
        // }

        // calculate internal loads
        assemble(config, geometry, beam_system);

        // step central differences (Includes updating displacements and nodal triads/Quaternions)
        // if (!set_disp)
        // {
        step_central_differences(config.dt, beam_system.u, beam_system.v, beam_system.M_inv,
                                 beam_system.J_u, beam_system.R, beam_system.R_static);
        // }
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
    assert(dx_min > 0);                     // just a check
    Scalar c = sqrt(config.E / config.rho); /*Speed of sound*/
    Scalar CFL = config.CFL;
    assert(CFL < 1 && CFL > 0);
    config.dt = CFL * dx_min / c;
}