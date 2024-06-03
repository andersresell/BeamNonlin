#include "../include/Solver.hpp"
#include "../include/HoleContact.hpp"
#include "../include/Numerics.hpp"
#include "../include/UserFunction.hpp"

Index n_glob;
Scalar t_glob;

static void calc_dt(Config &config, const Geometry &geometry, const Borehole &borehole);

static void calc_static_loads(const Config &config, const Geometry &geometry, vector<Vec3> &R_static_trans,
                              vector<Vec3> &R_static_rot);
static void check_energy_balance(const Config &config, const BeamSystem &beam_sys);
static void set_initial_configuration(const Config &config, vector<Vec3> &X, vector<Vec3> &d_trans,
                                      vector<Quaternion> &d_rot);

void solve(Config &config, Geometry &geometry, const Borehole &borehole) {

    create_output_dir(config);

    write_borehole(config, borehole);

    BeamSystem beam{config, geometry};

    set_initial_configuration(config, geometry.get_X(), beam.d_trans, beam.d_rot);

    calc_dt(config, geometry, borehole);

    initialize_hole_contact(config, geometry, borehole, beam);

    const Index n_steps = config.get_n_steps();

    calc_static_loads(config, geometry, beam.R_static_trans, beam.R_static_rot);

    calc_initial_accelerations(config, geometry, borehole, beam);

    printf("\n"
           "-----------------Starting simulation----------------\n"
           "The problem has Ne = %i elements\n"
           "Running with dt = %f for %i timesteps,\n"
           "for a total time of %f seconds\n"
           "--------------------------------------------------------\n",
           geometry.get_Ne(), config.dt, n_steps, config.dt * n_steps);

    Timer timer;

    timer.start_counter();

    for (Index n = 0; n <= n_steps; n++) {

        config.t = n * config.dt;
        config.n = n;

        // remove later
        n_glob = n;
        t_glob = config.t;
        DEBUG_ONLY(printf("n = %i", n));

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
        geometry.calc_cross_section_properties(ie, A, I_2, I_3,
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

static void check_energy_balance(const Config &config, const BeamSystem &beam) {
    if (!config.check_energy_balance) {
        return;
    }
    const Scalar KE = beam.KE;
    const Scalar W_int = beam.W_int;
    const Scalar W_ext = beam.W_ext;
    const Scalar tol = config.energy_balance_tol;
    const Scalar E_residual = KE + W_int - W_ext;

    if (abs(E_residual) > tol * max(KE, max(W_int, W_ext))) {
        printf("Warning: Energy balance is not obeyed, energy residual = %f\n", E_residual);
    } else if (isnan(E_residual)) {
        printf("Nan detected in energy residual\n");
    }
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
