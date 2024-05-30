#include "../include/Solver.hpp"
#include "../include/Numerics.hpp"

Index n_glob;
Scalar t_glob;

static void set_initial_configuration(const Config &config, vector<Vec3> &X, vector<Vec3> &d_trans,
                                      vector<Quaternion> &d_rot);
static void calc_dt(Config &config, const Geometry &geometry, const Borehole &borehole);

static void calc_static_loads(const Config &config, const Geometry &geometry, vector<Vec3> &R_static_trans,
                              vector<Vec3> &R_static_rot);
static void check_energy_balance(const Config &config, const BeamSystem &beam_sys);

void solve(Config &config, Geometry &geometry, const Borehole &borehole) {

    create_output_dir(config);

    write_borehole(config, borehole);

    BeamSystem beam{config, geometry};

    set_initial_configuration(config, geometry.get_X(), beam.d_trans, beam.d_rot);

    calc_dt(config, geometry, borehole);

    initialize_hole_contact(config, geometry, borehole, beam);

    const Index n_steps = config.get_n_steps();

    calc_static_loads(config, geometry, beam.R_static_trans, beam.R_static_rot);

    printf("\n"
           "-----------------Starting simulation----------------\n"
           "Running with dt = %f for %i timesteps,\n"
           "for a total time of %f seconds\n"
           "--------------------------------------------------------\n",
           config.dt, n_steps, config.dt * n_steps);

    Timer timer;

    timer.start_counter();

    for (Index n = 0; n <= n_steps; n++) {

        config.t = n * config.dt;
        config.n = n;

        // remove later
        n_glob = n;
        t_glob = config.t;

        if (n % 1000 == 0) {
            printf("\n---------------------------------\n"
                   "n = %i, t = %.5f"
                   "\n---------------------------------\n",
                   n, config.t);

            check_energy_balance(config, beam);
        }

        // save output
        save_csv(config, geometry, beam);

        // Debug stuff
        // constexpr bool set_disp = false;
        // if (set_disp && n == 0) {
        //     Index N = geometry.get_N();

        //     assert(geometry.get_N() == 2);
        //     // beam.u[1].trans = {0, u_new, 0};

        //     Scalar theta = 1 * M_PI / 180; // 1 deg
        //     Scalar c = cos(theta);
        //     Scalar s = sin(theta);

        //     // set external forces to zero
        //     for (Index i = 0; i < N; i++) {
        //         beam.R_ext_rot[i].setZero();
        //         beam.R_ext_trans[i].setZero();
        //         beam.R_static_rot[i].setZero();
        //         beam.R_static_trans[i].setZero();
        //     }

        //     Mat3 U = Mat3::Identity();
        //     Mat3 T = U;
        //     T = triad_from_euler_angles(10 * M_PI / 180, 0, 0);
        //     T = triad_from_euler_angles(0, 0, 45 * M_PI / 180);

        //     beam.d_trans[0] = Vec3{0, 0, 0};
        //     beam.d_trans[1] = Vec3{0, 0, 0};

        //     // beam_sys.d_rot[0].from_matrix(U);
        //     // beam_sys.d_rot[N - 1].from_matrix(T);
        //     // beam_sys.v_rot[0] = Vec3::Zero();
        //     // cout << "omega_mag_orig " << omega_u.norm() << endl;
        //     beam.v_rot[N - 1] = {10, 100, 0};
        //     // beam_sys.v_trans[N - 1] = {1000, 0, 0};
        // }

        step_explicit(config, geometry, borehole, beam);
    }
    timer.stop_counter();
    timer.print_elapsed_time();
}

static void calc_dt(Config &config, const Geometry &geometry, const Borehole &borehole) {

    const Scalar c = sqrt(config.E / config.rho); /*Speed of sound*/
    Scalar dx_min = std::numeric_limits<Scalar>::max();
    Scalar dt_min = std::numeric_limits<Scalar>::max();

    for (Index ie = 0; ie < geometry.get_Ne(); ie++) {

        Scalar dx = geometry.dx_e(ie);

        Scalar A{0}, I_2{0}, I_3{0}, unused{0};
        geometry.get_cross_section_properties(ie, A, I_2, I_3,
                                              unused); // should check if this works for rectangular cross sections

        Scalar I = max(I_2, I_3);
        Scalar r_g = sqrt(I / A); /*Radius of gyration*/

        Scalar crit_1 = sqrt(3) * dx_min * dx_min / (12 * c * r_g);
        Scalar crit_2 = dx / c;
        dt_min = min(crit_1, crit_2);
        dx_min = min(dx_min, dx);
        /*Testing the stability crit from table 6.1 in Belytscho*/
    }
    /*Find the smallest dx of the borehole*/
    if (config.contact_enabled) {
        Scalar dx_min_borehole = std::numeric_limits<Scalar>::max();
        for (Index ie_bh = 0; ie_bh < borehole.get_N_hole_elements(); ie_bh++) {
            dx_min_borehole = min(dx_min_borehole, borehole.dx(ie_bh));
        }
        if (dx_min > dx_min_borehole) {
            throw runtime_error("The minimum drillstring dx (" + to_string(dx_min) + " m)" +
                                "must be smaller than the minumum borehole dx (" + to_string(dx_min_borehole) +
                                " m). " + "This is to ensure robustness of the contact algorithm\n");
        }
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

static void calc_static_loads(const Config &config, const Geometry &geometry, vector<Vec3> &R_static_trans,
                              vector<Vec3> &R_static_rot) {
    const Index N = geometry.get_N();
    const Index Ne = N - 1;

    /*Set to zero first*/
    assert(R_static_trans.size() == N && R_static_rot.size() == N);
    for (Index i = 0; i < N; i++) {
        R_static_trans[i].setZero();
        R_static_rot[i].setZero();
    }

    const bool gravity_enabled = config.gravity_enabled;
    const Vec3 &g = config.gravity_acc;
    const Scalar rho = config.rho;

    for (Index ie = 0; ie < Ne; ie++) {
        if (gravity_enabled) {
            const Scalar m = geometry.dx_e(ie) * geometry.A_e(ie) * rho;

            R_static_trans[ie] += 0.5 * m * g;
            R_static_trans[ie + 1] += 0.5 * m * g;
            // R_static[ie].trans.y() += m * STANDARD_GRAVITY / 2;
            // R_static[ie + 1].trans.y() += m * STANDARD_GRAVITY / 2;
        }
    }

    /*Point loads*/
    /*If specified the point loads will be applied relative to the fixed orientation of the first node (typically
     * cantilever case)*/
    Mat3 R_base =
        config.point_loads_rel_to_base_orientation ? config.bc_orientation_base.to_matrix() : Mat3::Identity();
    for (const PointLoad &pl : config.R_point_static) {
        R_static_trans[pl.i] += R_base * pl.load_trans;

        R_static_rot[pl.i] += R_base * pl.load_rot; // * sin(config.t * 10);
    }
}

static void check_energy_balance(const Config &config, const BeamSystem &beam_sys) {
    if (!config.check_energy_balance) {
        return;
    }
    const Scalar KE = beam_sys.KE;
    const Scalar W_int = beam_sys.W_int;
    const Scalar W_ext = beam_sys.W_ext;
    const Scalar tol = config.energy_balance_tol;
    const Scalar E_residal = KE + W_int - W_ext;

    if (abs(E_residal) > tol * max(KE, max(W_int, W_ext))) {
        printf("Warning: Energy balance is not obeyed, energy residual = %f\n", E_residal);
    } else if (isnan(E_residal)) {
        printf("Nan detected in energy residual\n");
    }
}

static void calc_element_forces_local_rotated_TEST(Scalar ri, Scalar ro, Scalar l0, Scalar E, Scalar G, Scalar ul,
                                                   Scalar theta_1l, Scalar theta_2l, Scalar theta_3l, Scalar theta_4l,
                                                   Scalar theta_5l, Scalar theta_6l, Vec3 &f1, Vec3 &m1, Vec3 &f2,
                                                   Vec3 &m2) {
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

static void set_initial_configuration(const Config &config, vector<Vec3> &X, vector<Vec3> &d_trans,
                                      vector<Quaternion> &d_rot) {

    assert(X.size() == d_trans.size() && X.size() == d_rot.size());
    if (config.bc_case == BC_Case::CANTILEVER) {
        /*Rotate the reference configuration rigidly so that it matches The orientation at the base. Also set all
        rotations equal to the base rotation*/
        const Quaternion &q_base = config.bc_orientation_base;
        const Mat3 R = q_base.to_matrix();
        assert(X[0].isApprox(Vec3::Zero()));
        for (Index i = 0; i < X.size(); i++) {
            X[i] = R * X[i];
            d_rot[i] = q_base;
        }
    }
}