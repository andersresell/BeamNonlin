#include "../include/Solver.hpp"

void solve(Config &config, const Geometry &geometry)
{

    create_output_dir(config);
    BeamSystem beam_sys{config, geometry};
    calc_dt(config, geometry);
    const Index n_steps = config.get_n_steps();

    calc_static_loads(config, geometry, beam_sys.R_static_trans, beam_sys.R_static_rot);

    printf("\n"
           "-----------------Starting simulation----------------\n"
           "Running with dt = %f for %i timesteps,\n"
           "for a total time of %f seconds\n"
           "--------------------------------------------------------\n",
           config.dt, n_steps, config.dt * n_steps);

    assemble(config, geometry, beam_sys); //?

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
            printf("\n---------------------------------\n"
                   "n = %i, t = %.5f"
                   "\n---------------------------------\n",
                   n, config.t);

            check_energy_balance(config, beam_sys);
        }

        // save output
        save_csv(config, geometry, beam_sys);

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
        step_explicit(config, geometry, beam_sys);

        // calculate internal loads
        // assemble(config, geometry, beam_system);

        // step central differences (Includes updating displacements and nodal triads/Quaternions)
        // if (!set_disp)
        // // {
        // step_central_differences(config.dt, N, beam_system.u.data(), beam_system.v.data(), beam_system.M.data(),
        //                          beam_system.J_u.data(), beam_system.R_int.data(), beam_system.R_ext.data(),
        //                          check_energy_balance, config.W_int, config.W_ext, config.KE);
        // }
        // set_simple_bc(config.bc_case, geometry, beam_system);
    }
    timer.stop_counter();
    timer.print_elapsed_time();
}

void calc_dt(Config &config, const Geometry &geometry)
{
    Scalar c = sqrt(config.E / config.rho); /*Speed of sound*/
    Scalar dx_min = std::numeric_limits<Scalar>::max();
    Scalar dt_min = std::numeric_limits<Scalar>::max();
    for (Index ie = 0; ie < geometry.get_Ne(); ie++)
    {
        Scalar dx = geometry.dx_e(ie);
        Scalar I = geometry.I_e(ie);
        Scalar A = geometry.A_e(ie);
        Scalar r_g = sqrt(I / A); /*Radius of gyration*/

        Scalar crit_1 = sqrt(3) * dx_min * dx_min / (12 * c * r_g);
        Scalar crit_2 = dx / c;
        dt_min = min(crit_1, crit_2);
        dx_min = min(dx_min, dx);
        /*Testing the stability crit from table 6.1 in Belytscho*/
    }
    assert(dx_min > 0); // just a check
    assert(config.CFL < 1 && config.CFL > 0);
    config.dt = config.CFL * dt_min;
    cout << "----------------------- Choosing dt ----------------------\n"
         << "dt from table 6.1 Belytcho: " << dt_min << endl
         << "dt from l_min/c: " << dx_min / c << endl
         << "dt is chosen as CFL * dt_min, where CFL = " << config.CFL << endl
         << "dt = " << config.dt << endl
         << "----------------------------------------------------------\n";
}

void calc_static_loads(const Config &config, const Geometry &geometry,
                       vector<Vec3> &R_static_trans, vector<Vec3> &R_static_rot)
{
    const Index N = geometry.get_N();
    const Index Ne = N - 1;

    /*Set to zero first*/
    assert(R_static_trans.size() == N && R_static_rot.size() == N);
    for (Index i = 0; i < N; i++)
    {
        R_static_trans[i].setZero();
        R_static_rot[i].setZero();
    }

    const bool gravity_enabled = config.gravity_enabled;
    const Vec3 &g = config.gravity_acc;
    const Scalar rho = config.rho;

    for (Index ie = 0; ie < Ne; ie++)
    {
        if (gravity_enabled)
        {
            const Scalar m = geometry.dx_e(ie) * geometry.A_e(ie) * rho;
            R_static_trans[ie] += 0.5 * m * g;
            R_static_trans[ie + 1] += 0.5 * m * g;
            // R_static[ie].trans.y() += m * STANDARD_GRAVITY / 2;
            // R_static[ie + 1].trans.y() += m * STANDARD_GRAVITY / 2;
        }
    }

    /*Point loads*/
    for (const PointLoad &pl : config.R_point_static)
    {
        R_static_trans[pl.i] += pl.load_trans;
        R_static_rot[pl.i] += pl.load_rot;
    }
}