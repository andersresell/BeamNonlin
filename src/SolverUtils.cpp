#include "../include/SolverUtils.hpp"

BeamSystem::BeamSystem(const Config &config, const Geometry &geometry) : W_int{0}, W_ext{0}, KE{0} {
    const Index N = geometry.get_N();
    const Index Ne = N - 1;
    d_trans.resize(N, Vec3::Zero());
    d_rot.resize(N);
    v_trans.resize(N, Vec3::Zero());
    v_rot.resize(N, Vec3::Zero());
    a_rot.resize(N, Vec3::Zero());
    a_trans.resize(N, Vec3::Zero());
    L_rot.resize(N, Vec3::Zero());
    m_rot.resize(N, Vec3::Zero());
    R_int_trans.resize(N, Vec3::Zero());
    R_int_rot.resize(N, Vec3::Zero());
    R_ext_trans.resize(N, Vec3::Zero());
    R_ext_rot.resize(N, Vec3::Zero());
    R_static_trans.resize(N, Vec3::Zero());
    R_static_rot.resize(N, Vec3::Zero());
    M.resize(N, 0.0);
    J_u.resize(N, Vec3::Zero());

    if (config.contact_enabled) {
        hole_index.resize(N);
    }

    if (config.check_energy_balance) {
        delta_d_trans.resize(N);
        delta_d_rot.resize(N);
    }

    for (Index i = 0; i < N; i++) {
        d_rot[i].from_matrix(Mat3::Identity());
    }

    const Scalar rho = config.rho;

    for (Index ie = 0; ie < Ne; ie++) {
        const Scalar dx = geometry.dx_e(ie);
        const Scalar A = geometry.A_e(ie);
        const Scalar m = rho * dx * A;

        /*Lumped mass*/

        M[ie] += m / 2;
        M[ie + 1] += m / 2;

        /*Moment of inertia*/

        /*The reference configuration of the beam is along the x axis, with
        the local element frame E = identity*/

        /*Moment of inertia tensor of a hollow cylinder that lies along
        the x axis.

        Ju = [[J11 0    0 ]
              [0   J22  0 ]
              [0   0   J33]]

        J11 = integral rho (y² + z²)dV
        J22 = integral rho (x² + z²)dV
        J44 = integral rho (x² + y²)dV

        from https://en.wikipedia.org/wiki/List_of_moments_of_inertia :
        */
        Scalar Je11{};
        switch (geometry.get_cross_section_type()) {
        case CrossSectiontype::PIPE: {
            const Scalar ro = geometry.ro_e(ie);
            const Scalar ri = geometry.ri_e(ie);
            Je11 = rho * M_PI * dx / 2 * (powi<4>(ro) - powi<4>(ri));
            break;
        }
        case CrossSectiontype::RECANGLE: {
            const Scalar h2 = geometry.h2_e(ie);
            const Scalar h3 = geometry.h3_e(ie);
            Je11 = m / 12 * (h2 * h2 + h3 * h3);
            break;
        }
        default:
            assert(false);
        }
        // Je11 *= 30;

        const Scalar Je22 = 1.0 / 12 * m * dx * dx; // Moment of inertia of thin rod. Using this instead of the exact
                                                    // moment of inertia as Belytscho does Je22 *= 1;
        const Scalar Je33 = Je22;

        // Not the exact same procedure as proposed in Crisfield.
        // Seems like the moment of inertia is 1/12*m*l² there for y and z, i.e a thin rod

        /*Note: For now I assume that the reference configuration is a single straight beam oriented along
        the global x axis. If this is changed in the future I will need to implement the procedure in
        Crisfield given by eqns 24.141 and 24.142 which involves an eigenvalue decomposition to find
        the diagonal inertia matrix for each node. Since the reference configuration used now already
        leads to diagonal system, this wont make any difference.*/
        const Vec3 Je = {Je11, Je22, Je33};
        J_u[ie] += 0.5 * Je;
        J_u[ie + 1] += 0.5 * Je;
    }

    // /*Inverting lumped mass*/
    // for (Index i = 0; i < N; i++)
    // {
    //     M_inv[i] = 1 / M_inv[i];
    // }
}

void save_csv(const Config &config, const Geometry &geometry, const BeamSystem &beam) {
    using namespace std;
    const Index n = config.n;
    const Index n_w = config.n_write;
    if (!config.save_csv || n % n_w != 0)
        return;

    string filename = config.output_dir + to_string(n) + ".csv";
    ofstream ost{filename};

    if (!ost) {
        throw runtime_error{"Failed to open output file: " + filename + "\n"};
    }

    const Index N = geometry.get_N();

    /*Create "header"*/
    ost << "#N, n_steps, n_write, t, dt, check_energy_balance, contact_enabled\n"
        << N << "," << config.get_n_steps() << "," << n_w << "," << config.t << "," << config.dt << ","
        << config.check_energy_balance << "," << config.contact_enabled << "\n";
    /*Write energy*/
    if (config.check_energy_balance) {
        const Scalar E_res = beam.KE + beam.W_int - beam.W_ext;
        ost << "KE, W_int, W_ext, E_res\n"
            << beam.KE << ", " << beam.W_int << ", " << beam.W_ext << ", " << E_res << "\n ";
    } else {
        ost << "\n\n";
    }

    /*Write solution*/
    Index w = 12;
    ost << setw(w) << "#X1," << setw(w) << "X2," << setw(w) << "X3," << setw(w) << "d1," << setw(w) << "d2," << setw(w)
        << "d3," << setw(w) << "U11," << setw(w) << "U12," << setw(w) << "U13," << setw(w) << "U21," << setw(w)
        << "U22," << setw(w) << "U23," << setw(w) << "U31," << setw(w) << "U32," << setw(w) << "U33 " << setw(w)
        << "v1," << setw(w) << "v2," << setw(w) << "v3," << setw(w) << "omega_u_1," << setw(w) << "omega_u_2,"
        << setw(w) << "omega_u_3"
        << "\n";
    w -= 1;
    for (Index i = 0; i < N; i++) {
        const Vec3 &X = geometry.get_X()[i];
        const Vec3 &d = beam.d_trans[i];
        const Mat3 &U = beam.d_rot[i].to_matrix();
        const Vec3 &v = beam.v_trans[i];
        const Vec3 &omega_u = beam.v_rot[i];

        ost << setw(w) << X.x() << "," << setw(w) << X.y() << "," << setw(w) << X.z() << "," << setw(w) << d.x() << ","
            << setw(w) << d.y() << "," << setw(w) << d.z() << "," << setw(w) << U(0, 0) << "," << setw(w) << U(0, 1)
            << "," << setw(w) << U(0, 2) << "," << setw(w) << U(1, 0) << "," << setw(w) << U(1, 1) << "," << setw(w)
            << U(1, 2) << "," << setw(w) << U(2, 0) << "," << setw(w) << U(2, 1) << "," << setw(w) << U(2, 2) << ","
            << setw(w) << v.x() << "," << setw(w) << v.y() << "," << setw(w) << v.z() << "," << setw(w) << omega_u.x()
            << "," << setw(w) << omega_u.y() << "," << setw(w) << omega_u.z() << "\n";
    }
}

void write_borehole(const Config &config, const Borehole &borehole) {
    using namespace std;
    if (!config.contact_enabled)
        return;

    string filename = config.output_dir + "borehole.csv";
    ofstream ost{filename};

    if (!ost) {
        throw runtime_error{"Failed to open file: " + filename + "\n"};
    }

    const Index N_hole = borehole.get_N_hole_nodes();

    /*Create "header"*/
    ost << "#N_hole \n" << N_hole << "\n";

    /*Write solution*/
    Index w = 12;
    ost << setw(w) << "#X1," << setw(w) << "X2," << setw(w) << "X3," << setw(w) << "r_hole"
        << "\n";
    w -= 1;
    for (Index i = 0; i < N_hole; i++) {
        const Vec3 &x = borehole.get_x()[i];
        constexpr Scalar NaN = std::numeric_limits<Scalar>::quiet_NaN(); // Setting the last undefined value to nan
                                                                         // (since there is one less r than x)
        const Scalar r = (i < N_hole - 1) ? borehole.get_r_hole_element()[i] : NaN;

        ost << setw(w) << x.x() << "," << setw(w) << x.y() << "," << setw(w) << x.z() << "," << setw(w) << r << "\n";
    }
}

void create_output_dir(Config &config) {
    using namespace std;

    string &base_dir = config.base_dir;
    string &output_dir = config.output_dir;

    assert(base_dir != "");
    if (base_dir[base_dir.size() - 1] != '/') /*Adds a end slash if it is missing from the base_dir path*/
        base_dir += '/';

    output_dir = base_dir + "output/";

    if (filesystem::exists(output_dir)) {
        if (!filesystem::remove_all(output_dir)) {
            throw runtime_error{"Failed to remove old output directory: " + output_dir};
        }
    }
    if (!filesystem::create_directory(output_dir)) {
        throw runtime_error{string("Failed to create output directory:" + output_dir + "\n")};
    }
}
