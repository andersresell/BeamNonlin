#include "../include/SolverUtils.hpp"

BeamSystem::BeamSystem(const Config &config, const Geometry &geometry)
{
    Index N = geometry.get_N();
    Index Ne = N - 1;
    u.resize(N);
    v.resize(N);
    R.resize(N);
    R_static.resize(N);
    M_inv.resize(N);
    J_u.resize(N, Vec3::Zero());

    const Scalar rho = config.rho;

    for (Index ie = 0; ie < Ne; ie++)
    {
        Scalar dx = geometry.dx_e(ie);
        Scalar A = geometry.A_e(ie);
        Scalar m = rho * dx * A;
        Scalar ro = geometry.ro_e(ie);
        Scalar ri = geometry.ri_e(ie);

        /*Lumped mass*/

        M_inv[ie] += m / 2;
        M_inv[ie + 1] += m / 2;

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
        Scalar Je11 = rho * M_PI * dx / 2 * (pow(ro, 4) - pow(ri, 4));
        Scalar Je22;
        // Je22 = rho * M_PI * dx / 12 * (3 * (pow(ro, 4) - pow(ri, 4)) + dx * dx * (ro * ro - ri * ri));
        Je22 = 1.0 / 12 * m * dx * dx; // Moment of inertia of thin rod. Using this instead of the exact moment of inertia as Belytscho does
        // Je22 *= 1;
        Scalar Je33 = Je22;

        // Not the exact same procedure as proposed in Crisfield.
        // Seems like the moment of inertia is 1/12*m*l² there for y and z, i.e a thin rod

        /*Note: For now I assume that the reference configuration is a single straight beam oriented along
        the global x axis. If this is changed in the future I will need to implement the procedure in
        Crisfield given by eqns 24.141 and 24.142 which involves an eigenvalue decomposition to find
        the diagonal inertia matrix for each node. Since the reference configuration used now already
        leads to diagonal system, this wont make any difference.*/
        Vec3 Je = {Je11, Je22, Je33};
        J_u[ie] += 0.5 * Je;
        J_u[ie + 1] += 0.5 * Je;
    }

    /*Inverting lumped mass*/
    for (Index i = 0; i < N; i++)
    {
        M_inv[i] = 1 / M_inv[i];
    }
}

void save_csv(const Config &config, const Geometry &geometry, const BeamSystem &beam_system)
{
    using namespace std;
    Index n = config.n;
    Index n_w = config.n_write;
    if (!config.save_csv || n % n_w != 0)
        return;

    string filename = config.output_dir + to_string(n) + ".csv";
    ofstream ost{filename};

    if (!ost)
    {
        throw runtime_error{"Failed to open output file: " + filename + "\n"};
    }

    Index N = geometry.get_N();
    Index n_steps = config.get_n_steps();
    Scalar t = config.t;
    Scalar dt = config.dt;

    /*Create "header"*/
    ost
        << "#N, n_steps, n_write, t, dt\n"
        << N << "," << n_steps << "," << n_w << "," << t << "," << dt << "\n"
        << "\n";
    /*Write solution*/
    Index w = 12;
    ost << setw(w) << "#X1," << setw(w) << "X2," << setw(w) << "X3,"
        << setw(w) << "u1," << setw(w) << "u2," << setw(w) << "u3,"
        << setw(w) << "U11," << setw(w) << "U12," << setw(w) << "U13,"
        << setw(w) << "U21," << setw(w) << "U22," << setw(w) << "U23,"
        << setw(w) << "U31," << setw(w) << "U32," << setw(w) << "U33 "
        << setw(w) << "v1," << setw(w) << "v2," << setw(w) << "v3,"
        << setw(w) << "omega_u_1," << setw(w) << "omega_u_2," << setw(w) << "omega_u_3"
        << "\n";
    w -= 1;
    for (Index i = 0; i < N; i++)
    {
        const Vec3 &X = geometry.get_X()[i];
        const Vec3 &u = beam_system.u[i].trans;
        const Mat3 &U = beam_system.u[i].rot.to_matrix();
        const Vec3 &v = beam_system.v[i].trans;
        const Vec3 &omega_u = beam_system.v[i].rot;

        ost << setw(w) << X.x() << "," << setw(w) << X.y() << "," << setw(w) << X.z() << ","
            << setw(w) << u.x() << "," << setw(w) << u.y() << "," << setw(w) << u.z() << ","
            << setw(w) << U(0, 0) << "," << setw(w) << U(0, 1) << "," << setw(w) << U(0, 2) << ","
            << setw(w) << U(1, 0) << "," << setw(w) << U(1, 1) << "," << setw(w) << U(1, 2) << ","
            << setw(w) << U(2, 0) << "," << setw(w) << U(2, 1) << "," << setw(w) << U(2, 2) << ","
            << setw(w) << v.x() << "," << setw(w) << u.y() << "," << setw(w) << u.z() << ","
            << setw(w) << omega_u.x() << "," << setw(w) << omega_u.y() << "," << setw(w) << omega_u.z() << "\n";
    }
}

void create_output_dir(Config &config)
{
    using namespace std;

    string &base_dir = config.base_dir;
    string &output_dir = config.output_dir;

    assert(base_dir != "");
    if (base_dir[base_dir.size() - 1] != '/') /*Adds a end slash if it is missing from the base_dir path*/
        base_dir += '/';

    output_dir = base_dir + "output/";

    if (filesystem::exists(output_dir))
    {
        if (!filesystem::remove_all(output_dir))
        {
            throw runtime_error{"Couldn't remove old output directory: " + output_dir};
        }
    }
    if (!filesystem::create_directory(output_dir))
    {
        throw runtime_error{string("Failed to create output directory:" + output_dir + "\n")};
    }
}
