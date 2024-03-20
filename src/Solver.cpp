#include "../include/Solver.hpp"
#include "Solver.inl"

void solve(Config &config, Geometry &geometry, const Borehole &borehole)
{

    create_output_dir(config);
    BeamSystem beam_sys{config, geometry};

    set_initial_configuration(config, geometry.get_X(), beam_sys.d_trans, beam_sys.d_rot);

    calc_dt(config, geometry);
    const Index n_steps = config.get_n_steps();

    calc_static_loads(config, geometry, beam_sys.R_static_trans, beam_sys.R_static_rot);

    printf("\n"
           "-----------------Starting simulation----------------\n"
           "Running with dt = %f for %i timesteps,\n"
           "for a total time of %f seconds\n"
           "--------------------------------------------------------\n",
           config.dt, n_steps, config.dt * n_steps);

    Timer timer;

    timer.start_counter();

    for (Index n = 0; n <= n_steps; n++)
    {

        config.t = n * config.dt;
        config.n = n;

        calc_static_loads(config, geometry, beam_sys.R_static_trans, beam_sys.R_static_rot);

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

        // Debug stuff
        constexpr bool set_disp = false;
        if (set_disp && n == 0)
        {
            Index N = geometry.get_N();

            assert(geometry.get_N() == 2);
            // beam_system.u[1].trans = {0, u_new, 0};

            Scalar theta = 1 * M_PI / 180; // 1 deg
            Scalar c = cos(theta);
            Scalar s = sin(theta);

            // set external forces to zero
            for (Index i = 0; i < N; i++)
            {
                beam_sys.R_ext_rot[i].setZero();
                beam_sys.R_ext_trans[i].setZero();
                beam_sys.R_static_rot[i].setZero();
                beam_sys.R_static_trans[i].setZero();
            }

            Mat3 U = Mat3::Identity();
            Mat3 T = U;
            T = triad_from_euler_angles(10 * M_PI / 180, 0, 0);
            T = triad_from_euler_angles(0, 0, 45 * M_PI / 180);

            beam_sys.d_trans[0] = Vec3{0, 0, 0};
            beam_sys.d_trans[1] = Vec3{0, 0, 0};

            // beam_sys.d_rot[0].from_matrix(U);
            // beam_sys.d_rot[N - 1].from_matrix(T);
            // beam_sys.v_rot[0] = Vec3::Zero();
            // cout << "omega_mag_orig " << omega_u.norm() << endl;
            beam_sys.v_rot[N - 1] = {10000, 0, 0};
            // beam_sys.v_trans[N - 1] = {1000, 0, 0};
        }

        step_explicit_NMB(config, geometry, borehole, beam_sys);

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

void step_explicit_NMB(Config &config, const Geometry &geometry, const Borehole &borehole, BeamSystem &beam_sys)
{
    const Scalar dt = config.dt;
    const Index N = geometry.get_N();
    vector<Vec3> &d_trans = beam_sys.d_trans;
    vector<Quaternion> &d_rot = beam_sys.d_rot;
    vector<Vec3> &v_trans = beam_sys.v_trans;
    vector<Vec3> &v_rot = beam_sys.v_rot;
    vector<Vec3> &a_trans = beam_sys.a_trans;
    vector<Vec3> &a_rot = beam_sys.a_rot;
    vector<Vec3> &L_rot = beam_sys.L_rot;
    vector<Vec3> &m_rot = beam_sys.m_rot;
    vector<Vec3> &R_int_trans = beam_sys.R_int_trans;
    vector<Vec3> &R_int_rot = beam_sys.R_int_rot;
    vector<Vec3> &R_ext_trans = beam_sys.R_ext_trans;
    vector<Vec3> &R_ext_rot = beam_sys.R_ext_rot;
    const vector<Scalar> &M = beam_sys.M;
    const vector<Vec3> &J_u = beam_sys.J_u;
    const bool check_energy_balance = config.check_energy_balance;
    const bool rayleigh_damping_enabled = config.rayleigh_damping_enabled;
    Scalar &W_int = beam_sys.W_int;
    Scalar &W_ext = beam_sys.W_ext;
    Scalar &KE = beam_sys.KE;
    vector<Vec3> &delta_d_trans = beam_sys.delta_d_trans; /*Only used if energy balance is checked*/
    vector<Vec3> &delta_d_rot = beam_sys.delta_d_rot;     /*Only used if energy balance is checked*/
    if (check_energy_balance)
    {
        assert(delta_d_trans.size() == N && delta_d_rot.size() == N);
    }
    else
    {
        assert(delta_d_trans.size() == 0 && delta_d_rot.size() == 0);
    }

    for (Index i = 0; i < N; i++)
    {
        const Vec3 delta_d = dt * v_trans[i] + 0.5 * dt * dt * a_trans[i];
        d_trans[i] += delta_d;
        if (check_energy_balance)
        {
            delta_d_trans[i] = delta_d;
        }
    }

    // rotations:
    for (Index i = 0; i < N; i++)
    {
        Quaternion &q = d_rot[i];
        const Mat3 &J = J_u[i].asDiagonal();
        Vec3 &omega_u = v_rot[i];
        Vec3 &alpha_u = a_rot[i];
        L_rot[i] = q.rotate_vector(J * omega_u); // Storing angular momentum L = U*J_u*omega_u at t_n for velocity update
        const Vec3 theta_u = dt * omega_u + 0.5 * dt * dt * alpha_u;
        q.exponential_map_body_frame(theta_u); // Update the rotation as U_{n+1} = U_n * exp(S(theta_u))

        if (check_energy_balance)
        {
            delta_d_rot[i] = q.rotate_vector(theta_u); // delta_d is stored in inertial frame
        }
    }

    /*Enforcing boundary conditions*/
    set_simple_bc(config, geometry, beam_sys);

    if (check_energy_balance)
    {
        work_update_partial(N, delta_d_trans.data(), delta_d_rot.data(), R_int_trans.data(),
                            R_int_rot.data(), R_ext_trans.data(), R_ext_rot.data(), W_ext, W_int);
    }

    /*Update internal and external forces*/
    assemble(config, geometry, beam_sys);
    if (rayleigh_damping_enabled)
    {
        add_mass_proportional_rayleigh_damping(N, config.alpha_rayleigh, M.data(), v_trans.data(), R_int_trans.data(),
                                               J_u.data(), d_rot.data(), v_rot.data(), R_int_rot.data());
    }

    if (check_energy_balance)
    {
        work_update_partial(N, delta_d_trans.data(), delta_d_rot.data(), R_int_trans.data(),
                            R_int_rot.data(), R_ext_trans.data(), R_ext_rot.data(), W_ext, W_int);
    }

    /*Update translation velocities and the translational accelerations */
    for (Index i = 0; i < N; i++)
    {
        const Vec3 a_trans_new = (R_ext_trans[i] - R_int_trans[i]) / M[i];
        // v_trans[i] += 0.5 * dt * (a_trans[i] + a_trans_new);
        v_trans[i] += dt * 0.5 * (a_trans[i] + a_trans_new);
        a_trans[i] = a_trans_new;
    }

    // rotations: Newmark body
    for (Index i = 0; i < N; i++)
    {
        const Quaternion &q = d_rot[i];
        const Mat3 &J = J_u[i].asDiagonal();
        const Vec3 &L_n = L_rot[i];
        Vec3 &omega_u = v_rot[i];
        Vec3 &alpha_u = a_rot[i];

        const Mat3 U_np = q.to_matrix(); // maybe optimize later
        const Vec3 m = R_ext_rot[i] - R_int_rot[i];

        const Vec3 omega_u_old = omega_u;
        omega_u = J.inverse() * U_np.transpose() * (L_n + dt * m_half);

        alpha_u = 2 / dt * (omega_u - omega_u_old) - alpha_u; // beta = 0.5 (not recommended by Simo, but seems to work better)
        // alpha_u = (omega_u - omega_u_old) / dt; // Update angular acceration in body frame
    }

    /*Enforcing boundary conditions*/
    set_simple_bc(config, geometry, beam_sys);

    if (check_energy_balance)
    {
        kinetic_energy_update(N, M.data(), J_u.data(), v_trans.data(), v_rot.data(), KE);
    }
}

void step_explicit_SW(Config &config, const Geometry &geometry, const Borehole &borehole, BeamSystem &beam_sys)
{
    const Scalar dt = config.dt;
    const Index N = geometry.get_N();
    vector<Vec3> &d_trans = beam_sys.d_trans;
    vector<Quaternion> &d_rot = beam_sys.d_rot;
    vector<Vec3> &v_trans = beam_sys.v_trans;
    vector<Vec3> &v_rot = beam_sys.v_rot;
    vector<Vec3> &a_trans = beam_sys.a_trans;
    vector<Vec3> &a_rot = beam_sys.a_rot;
    vector<Vec3> &L_rot = beam_sys.L_rot;
    vector<Vec3> &m_rot = beam_sys.m_rot;
    vector<Vec3> &R_int_trans = beam_sys.R_int_trans;
    vector<Vec3> &R_int_rot = beam_sys.R_int_rot;
    vector<Vec3> &R_ext_trans = beam_sys.R_ext_trans;
    vector<Vec3> &R_ext_rot = beam_sys.R_ext_rot;
    const vector<Scalar> &M = beam_sys.M;
    const vector<Vec3> &J_u = beam_sys.J_u;
    const bool check_energy_balance = config.check_energy_balance;
    const bool rayleigh_damping_enabled = config.rayleigh_damping_enabled;
    Scalar &W_int = beam_sys.W_int;
    Scalar &W_ext = beam_sys.W_ext;
    Scalar &KE = beam_sys.KE;
    vector<Vec3> &delta_d_trans = beam_sys.delta_d_trans; /*Only used if energy balance is checked*/
    vector<Vec3> &delta_d_rot = beam_sys.delta_d_rot;     /*Only used if energy balance is checked*/
    if (check_energy_balance)
    {
        assert(delta_d_trans.size() == N && delta_d_rot.size() == N);
    }
    else
    {
        assert(delta_d_trans.size() == 0 && delta_d_rot.size() == 0);
    }

    // bool half_step = false;

    // if (half_step)

    //     for (Index i = 0; i < N; i++)
    //     {
    //         v_trans[i] += dt / M[i] * (R_ext_trans[i] - R_int_trans[i]);
    //         delta_d_trans[i] = dt * v_trans[i];
    //     }
    // else
    // newmark beta with beta=0, gamma=1/2 (central difference)

    for (Index i = 0; i < N; i++)
    {
        const Vec3 delta_d = dt * v_trans[i] + 0.5 * dt * dt * a_trans[i];
        d_trans[i] += delta_d;
        if (check_energy_balance)
        {
            delta_d_trans[i] = delta_d;
        }
    }

    // rotations: Simo and Wong algorithm
    for (Index i = 0; i < N; i++)
    {
        Quaternion &q = d_rot[i];
        const Mat3 &J = J_u[i].asDiagonal();
        Vec3 &omega_u = v_rot[i];
        Vec3 &alpha_u = a_rot[i];
        // const Mat3 U_n = q.to_matrix(); // Copy of orientation at t_n. Optimize this later maybe.
        // L_rot[i] = U_n * J * omega_u;   // Storing angular momentum at t_n for velocity update
        L_rot[i] = q.rotate_vector(J * omega_u); // Storing angular momentum L = U*J_u*omega_u at t_n for velocity update
        const Vec3 theta_u = dt * omega_u + 0.5 * dt * dt * alpha_u;
        q.exponential_map_body_frame(theta_u); // Update the rotation as U_{n+1} = U_n * exp(S(theta_u))

        if (check_energy_balance)
        {
            delta_d_rot[i] = q.rotate_vector(theta_u); // delta_d is stored in inertial frame
        }
    }

    /*Enforcing boundary conditions*/
    set_simple_bc(config, geometry, beam_sys);

    if (check_energy_balance)
    {
        work_update_partial(N, delta_d_trans.data(), delta_d_rot.data(), R_int_trans.data(),
                            R_int_rot.data(), R_ext_trans.data(), R_ext_rot.data(), W_ext, W_int);
    }

    /*Update internal and external forces*/
    assemble(config, geometry, beam_sys);
    if (rayleigh_damping_enabled)
    {
        add_mass_proportional_rayleigh_damping(N, config.alpha_rayleigh, M.data(), v_trans.data(), R_int_trans.data(),
                                               J_u.data(), d_rot.data(), v_rot.data(), R_int_rot.data());
    }

    if (check_energy_balance)
    {
        work_update_partial(N, delta_d_trans.data(), delta_d_rot.data(), R_int_trans.data(),
                            R_int_rot.data(), R_ext_trans.data(), R_ext_rot.data(), W_ext, W_int);
    }

    /*Update translation velocities and the translational accelerations */
    for (Index i = 0; i < N; i++)
    {
        const Vec3 a_trans_new = (R_ext_trans[i] - R_int_trans[i]) / M[i];
        // v_trans[i] += 0.5 * dt * (a_trans[i] + a_trans_new);
        v_trans[i] += dt * 0.5 * (a_trans[i] + a_trans_new);
        a_trans[i] = a_trans_new;
    }

    // rotations: Simo and Wong algorithm
    for (Index i = 0; i < N; i++)
    {
        const Quaternion &q = d_rot[i];
        const Mat3 &J = J_u[i].asDiagonal();
        const Vec3 &L_n = L_rot[i];
        Vec3 &omega_u = v_rot[i];
        Vec3 &alpha_u = a_rot[i];

        const Mat3 U_np = q.to_matrix(); // maybe optimize later
        Vec3 &m = m_rot[i];              // Moment at t_{n+1/2}
        const Vec3 m_np = R_ext_rot[i] - R_int_rot[i];
        /*Evaluate moment at t_{n+1/2} by trapezoidal rule, i.e m_{n+1/2} = 1/2*(m_{n} + m_{n+1}) */
        Vec3 m_half = 0.5 * (m + m_np);

        m = m_np; // Update moment
        // if (i == 1)
        //     cout << "m_half\n"
        //          << m_half << endl;

        const Vec3 omega_u_old = omega_u;
        omega_u = J.inverse() * U_np.transpose() * (L_n + dt * m_half);
        // cout << "omega_u\n " << omega_u << endl;
#ifndef NDEBUG
        Vec3 L_np = U_np * J * omega_u;
        Vec3 res = L_np - L_n - dt * m_half;
        assert(is_close(res.norm(), 0.0));
#endif
        alpha_u = 2 / dt * (omega_u - omega_u_old) - alpha_u; // beta = 0.5 (not recommended by Simo, but seems to work better)
        // alpha_u = (omega_u - omega_u_old) / dt; // Update angular acceration in body frame
    }

    /*Enforcing boundary conditions*/
    set_simple_bc(config, geometry, beam_sys);

    if (check_energy_balance)
    {
        kinetic_energy_update(N, M.data(), J_u.data(), v_trans.data(), v_rot.data(), KE);
    }
}

void step_explicit_old(Config &config, const Geometry &geometry, const Borehole &borehole, BeamSystem &beam_sys)
{
    const Scalar dt = config.dt;
    const Index N = geometry.get_N();
    vector<Vec3> &d_trans = beam_sys.d_trans;
    vector<Quaternion> &d_rot = beam_sys.d_rot;
    vector<Vec3> &v_trans = beam_sys.v_trans;
    vector<Vec3> &v_rot = beam_sys.v_rot;
    vector<Vec3> &R_int_trans = beam_sys.R_int_trans;
    vector<Vec3> &R_int_rot = beam_sys.R_int_rot;
    vector<Vec3> &R_ext_trans = beam_sys.R_ext_trans;
    vector<Vec3> &R_ext_rot = beam_sys.R_ext_rot;
    const vector<Scalar> &M = beam_sys.M;
    const vector<Vec3> &J_u = beam_sys.J_u;
    const bool check_energy_balance = config.check_energy_balance;
    const bool rayleigh_damping_enabled = config.rayleigh_damping_enabled;
    Scalar &W_int = beam_sys.W_int;
    Scalar &W_ext = beam_sys.W_ext;
    Scalar &KE = beam_sys.KE;
    vector<Vec3> &delta_d_trans = beam_sys.delta_d_trans; /*Only used if energy balance is checked*/
    vector<Vec3> &delta_d_rot = beam_sys.delta_d_rot;     /*Only used if energy balance is checked*/

    if (check_energy_balance)
    {
        assert(delta_d_trans.size() == N && delta_d_rot.size() == N);
    }
    else
    {
        assert(delta_d_trans.size() == 0 && delta_d_rot.size() == 0);
    }

    // /*velocity at t_{n+1/2}*/
    // velocity_update_partial(dt, N, M.data(), J_u.data(), R_int_trans.data(), R_int_rot.data(),
    //                         R_ext_trans.data(), R_ext_rot.data(), v_trans.data(), v_rot.data());

    /*Enforcing boundary conditions*/
    set_simple_bc(config, geometry, beam_sys);

    /*Checking energy balance in two steps.
    The equation to be computed is W_{n+1} = W_n + delta_u/2*(Rn + R_{n+1})
    This will be done in two steps in order to not need double storage for R
    First W += delta_u/2*R_n, then updating R, then W +=  * delta_u/2*R_{n+1}
    */
    if (check_energy_balance)
    {
        calc_delta_d(dt, N, delta_d_trans.data(), delta_d_rot.data(), v_trans.data(), v_rot.data());

        work_update_partial(N, delta_d_trans.data(), delta_d_rot.data(), R_int_trans.data(),
                            R_int_rot.data(), R_ext_trans.data(), R_ext_rot.data(), W_ext, W_int);
    }

    /*Displacement at t_{n+1}*/
    displacement_update(dt, N, v_trans.data(), v_rot.data(), d_trans.data(), d_rot.data());

    //////// Move this to a unfied force function
    /*Updating external and internal forces*/
    assemble(config, geometry, beam_sys);

    if (rayleigh_damping_enabled)
    {
        add_mass_proportional_rayleigh_damping(N, config.alpha_rayleigh, M.data(), v_trans.data(), R_int_trans.data(),
                                               J_u.data(), d_rot.data(), v_rot.data(), R_int_rot.data());
    }
    //////////

    /*Since the moments are allways used in the body frame, these are rotated
    once and for all instead of every time they are needed.*/
    rotate_moment_to_body_frame(N, d_rot.data(), R_int_rot.data(), R_ext_rot.data());

    if (check_energy_balance)
    {
        work_update_partial(N, delta_d_trans.data(), delta_d_rot.data(), R_int_trans.data(),
                            R_int_rot.data(), R_ext_trans.data(), R_ext_rot.data(), W_ext, W_int);
    }

    /*velocity at t_{n+1} */
    // velocity_update_partial(dt, N, M.data(), J_u.data(), R_int_trans.data(), R_int_rot.data(),
    //                         R_ext_trans.data(), R_ext_rot.data(), v_trans.data(), v_rot.data());

    /*Enforcing boundary conditions*/
    set_simple_bc(config, geometry, beam_sys);

    if (check_energy_balance)
    {
        kinetic_energy_update(N, M.data(), J_u.data(), v_trans.data(), v_rot.data(), KE);
    }
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
    /*If specified the point loads will be applied relative to the fixed orientation of the first node (typically cantilever case)*/
    Mat3 R_base = config.point_loads_rel_to_base_orientation ? config.bc_orientation_base.to_matrix() : Mat3::Identity();
    for (const PointLoad &pl : config.R_point_static)
    {
        R_static_trans[pl.i] += R_base * pl.load_trans;

        R_static_rot[pl.i] += R_base * pl.load_rot; // * sin(config.t * 10);
    }
}

void check_energy_balance(const Config &config, const BeamSystem &beam_sys)
{
    if (!config.check_energy_balance)
    {
        return;
    }
    const Scalar KE = beam_sys.KE;
    const Scalar W_int = beam_sys.W_int;
    const Scalar W_ext = beam_sys.W_ext;
    const Scalar tol = config.energy_balance_tol;
    const Scalar E_residal = KE + W_int - W_ext;

    if (abs(E_residal) > tol * max(KE, max(W_int, W_ext)))
    {
        printf("Warning: Energy balance is not obeyed, energy residual = %f\n", E_residal);
    }
    else if (isnan(E_residal))
    {
        printf("Nan detected in energy residual\n");
    }
}

inline void calc_element_forces_local_rotated_TEST(Scalar ri, Scalar ro, Scalar l0, Scalar E, Scalar G, Scalar ul,
                                                   Scalar theta_1l, Scalar theta_2l, Scalar theta_3l, Scalar theta_4l,
                                                   Scalar theta_5l, Scalar theta_6l, Vec3 &f1, Vec3 &m1, Vec3 &f2, Vec3 &m2)
{
    const Scalar A = M_PI * (ro * ro - ri * ri);
    const Scalar I = M_PI / 4 * (ro * ro * ro * ro - ri * ri * ri * ri);
    const Scalar J = 2 * I;

    /*Normal force (F1)
     [[F1], = A*E/l0[[ 1 -1],*[[0],
      [F4]]          [-1  1]]  [ul]]
    */
    const Scalar F1 = A * E * (-ul) / l0;
    const Scalar F4 = -F1;

    /*--------------------------------------------------------------------
    Torsion:
    Used the theory from this link:
    https://www.acs.psu.edu/drussell/Demos/Torsional/torsional.html
    which is a simple wave equation. As long as circular bars are used
    I_p = K and these terms disappear.
    Governing equation is then:
    I_p * rho * phitors_tt = G * K * phitors_xx
    --------------------------------------------------------------------*/
    const Scalar K = J; // cicular cross-section
    const Scalar M1 = G * K * (theta_1l - theta_4l) / l0;
    const Scalar M4 = -M1;

    /*Bending. Euler Bernoulli is used, symmetrical cross section
    The original matrix reads
    EI/L³ *[[ 12  6L  -12  6L ]
            [ 6L  4L² -6L  2L²]
            [-12 -6L   12 -6L ]
            [ 6L  2L² -6L  4L²]]
    multiplied by [w1 theta1, w2 theta2],
    however, it can be simplified, since w1=w2=0,
    rows 1 and 3 can be skipped resulting in
    EI/L²**[[ 6   6 ]*[[theta1 ]
            [ 4L  2L]  [theta2 ]]
            [-6  -6 ]
            [ 2L  4L]]
    */

    const Scalar F2 = E * I / (l0 * l0) * (6 * theta_2l + 6 * theta_5l);
    const Scalar M2 = E * I / l0 * (4 * theta_2l + 2 * theta_5l);
    const Scalar F5 = -F2;
    const Scalar M5 = E * I / l0 * (2 * theta_2l + 4 * theta_5l);

    const Scalar F3 = E * I / (l0 * l0) * (6 * theta_3l + 6 * theta_6l);
    const Scalar M3 = E * I / l0 * (4 * theta_3l + 2 * theta_6l);
    const Scalar F6 = -F3;
    const Scalar M6 = E * I / l0 * (2 * theta_3l + 4 * theta_6l);
    // const Scalar F2 = E * I / (l0 * l0) * (6 * theta_2l + 6 * theta_5l);
    // const Scalar M2 = E * I / (l0 * l0) * (4 * l0 * theta_2l + 2 * l0 * theta_5l);
    // const Scalar F5 = -F2;
    // const Scalar M5 = E * I / (l0 * l0) * (2 * l0 * theta_2l + 4 * l0 * theta_5l);
    f1.x() = F1;
    f1.y() = F2;
    f1.z() = F3;

    m1.x() = M1;
    m1.y() = -M3;
    m1.z() = M2;

    f2.x() = F4;
    f2.y() = F5;
    f2.z() = F6;

    m2.x() = M4;
    m2.y() = -M6;
    m2.z() = M5;
}

void set_initial_configuration(const Config &config, vector<Vec3> &X, vector<Vec3> &d_trans, vector<Quaternion> &d_rot)
{

    assert(X.size() == d_trans.size() && X.size() == d_rot.size());
    if (config.bc_case == BC_Case::CANTILEVER)
    {
        /*Rotate the reference configuration rigidly so that it matches The orientation at the base. Also set all
        rotations equal to the base rotation*/
        const Quaternion &q_base = config.bc_orientation_base;
        const Mat3 R = q_base.to_matrix();
        assert(X[0].isApprox(Vec3::Zero()));
        for (Index i = 0; i < X.size(); i++)
        {
            X[i] = R * X[i];
            d_rot[i] = q_base;
        }
    }
}