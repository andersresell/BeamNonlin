#include "../include/Numerics.hpp"
#include "../include/BattiniBeam.hpp"
#include "../include/Containers.hpp"
#include "../include/CrisfieldBeam.hpp"
#include "../include/HoleContact.hpp"
#include "../include/SolverUtils.hpp"
#include "../include/UserFunction.hpp"

static void calc_forces(const Config &config, const Geometry &geometry, const Borehole &borehole, BeamSystem &beam);

template <Index i_first>
static void calc_inner_forces(const Config &config, const Geometry &geometry, BeamSystem &beam);
static void set_simple_bc(const Config &config, const Geometry &geometry, BeamSystem &beam);
static void work_update_partial(Index N, const vector<Vec3> &delta_d_trans, const vector<Vec3> &delta_d_rot,
                                const vector<Vec3> &R_int_trans, const vector<Vec3> &R_int_rot,
                                const vector<Vec3> &R_ext_trans, const vector<Vec3> &R_ext_rot, Scalar &W_ext,
                                Scalar &W_int);
static void kinetic_energy_update(Index N, const vector<Scalar> &M, const vector<Vec3> &J_u,
                                  const vector<Vec3> &v_trans, const vector<Vec3> &v_rot, Scalar &KE);
static void zero_internal_and_set_static_forces(const Index N, const vector<Vec3> &R_static_trans,
                                                const vector<Vec3> &R_static_rot, vector<Vec3> &R_int_trans,
                                                vector<Vec3> &R_int_rot, vector<Vec3> &R_ext_trans,
                                                vector<Vec3> &R_ext_rot);
static void add_mass_proportional_rayleigh_damping(Index N, Scalar alpha, const vector<Scalar> &M,
                                                   const vector<Vec3> &v_trans, vector<Vec3> &R_int_trans,
                                                   const vector<Vec3> &J_u, const vector<Quaternion> &d_rot,
                                                   const vector<Vec3> &v_rot, vector<Vec3> &R_int_rot);

/*--------------------------------------------------------------------
Explicit solution by the Simo-Wong algorithm for rotations and
non-staggered central differences for translations
--------------------------------------------------------------------*/
void step_explicit(Config &config, const Geometry &geometry, const Borehole &borehole, BeamSystem &beam) {
    const Scalar dt = config.dt;
    const Index N = geometry.get_N();
    vector<Vec3> &d_trans = beam.d_trans;
    vector<Quaternion> &d_rot = beam.d_rot;
    vector<Vec3> &v_trans = beam.v_trans;
    vector<Vec3> &v_rot = beam.v_rot;
    vector<Vec3> &a_trans = beam.a_trans;
    vector<Vec3> &a_rot = beam.a_rot;
    vector<Vec3> &L_rot = beam.L_rot;
    vector<Vec3> &m_rot = beam.m_rot;
    vector<Vec3> &R_int_trans = beam.R_int_trans;
    vector<Vec3> &R_int_rot = beam.R_int_rot;
    vector<Vec3> &R_ext_trans = beam.R_ext_trans;
    vector<Vec3> &R_ext_rot = beam.R_ext_rot;
    const vector<Scalar> &M = beam.M;
    const vector<Vec3> &J_u = beam.J_u;
    const bool check_energy_balance = config.check_energy_balance;
    Scalar &W_int = beam.W_int;
    Scalar &W_ext = beam.W_ext;
    Scalar &KE = beam.KE;
    vector<Vec3> &delta_d_trans = beam.delta_d_trans; /*Only used if energy balance is checked*/
    vector<Vec3> &delta_d_rot = beam.delta_d_rot;     /*Only used if energy balance is checked*/
    if (check_energy_balance) {
        assert(delta_d_trans.size() == N && delta_d_rot.size() == N);
    } else {
        assert(delta_d_trans.size() == 0 && delta_d_rot.size() == 0);
    }

#pragma omp parallel for
    for (Index i = 0; i < N; i++) {
        const Vec3 delta_d = dt * v_trans[i] + 0.5 * dt * dt * a_trans[i];
        d_trans[i] += delta_d;
        if (check_energy_balance) {
            delta_d_trans[i] = delta_d;
        }
    }

    // rotations: Simo and Wong algorithm

#pragma omp parallel for
    for (Index i = 0; i < N; i++) {
        Quaternion &q = d_rot[i];
        const Mat3 &J = J_u[i].asDiagonal();
        Vec3 &omega_u = v_rot[i];
        Vec3 &alpha_u = a_rot[i];
        L_rot[i] =
            q.rotate_vector(J * omega_u); // Storing angular momentum L = U*J_u*omega_u at t_n for velocity update
        const Vec3 theta_u = dt * omega_u + 0.5 * dt * dt * alpha_u;
        q.exponential_map_body_frame(theta_u); // Update the rotation as U_{n+1} = U_n * exp(S(theta_u))

        if (check_energy_balance) {
            delta_d_rot[i] = q.rotate_vector(theta_u); // delta_d is stored in inertial frame
        }
    }

    /*Enforcing boundary conditions*/
    set_simple_bc(config, geometry, beam);

    if (check_energy_balance) {
        work_update_partial(N, delta_d_trans, delta_d_rot, R_int_trans, R_int_rot, R_ext_trans, R_ext_rot, W_ext,
                            W_int);
    }

    /*Update internal and external forces*/
    calc_forces(config, geometry, borehole, beam);

    if (check_energy_balance) {
        work_update_partial(N, delta_d_trans, delta_d_rot, R_int_trans, R_int_rot, R_ext_trans, R_ext_rot, W_ext,
                            W_int);
    }

    /*Update translation velocities and the translational accelerations */

#pragma omp parallel for
    for (Index i = 0; i < N; i++) {
        const Vec3 a_trans_new = (R_ext_trans[i] - R_int_trans[i]) / M[i];
        v_trans[i] += dt * 0.5 * (a_trans[i] + a_trans_new);
        a_trans[i] = a_trans_new;
    }

    // rotations: Simo and Wong algorithm

#pragma omp parallel for
    for (Index i = 0; i < N; i++) {
        const Mat3 &J = J_u[i].asDiagonal();
        const Vec3 &L_n = L_rot[i];
        Vec3 &omega_u = v_rot[i];
        Vec3 &alpha_u = a_rot[i];
        const Mat3 U_np = d_rot[i].to_matrix();
        Vec3 &m = m_rot[i]; // Moment at t_{n+1/2}
        const Vec3 m_np = R_ext_rot[i] - R_int_rot[i];
        /*Evaluate moment at t_{n+1/2} by trapezoidal rule, i.e m_{n+1/2} = 1/2*(m_{n} + m_{n+1}) */
        Vec3 m_half = 0.5 * (m + m_np);
        m = m_np; // Update moment
        const Vec3 omega_u_old = omega_u;
        omega_u = J.inverse() * U_np.transpose() * (L_n + dt * m_half);
#ifndef NDEBUG
        const Vec3 L_np = U_np * J * omega_u;
        const Vec3 res = L_np - L_n - dt * m_half;
        assert(is_close(res.norm(), 0.0));
#endif
        alpha_u = 2 / dt * (omega_u - omega_u_old) -
                  alpha_u; // beta = 0.5 (not recommended by Simo, but seems to work better)
        // alpha_u = (omega_u - omega_u_old) / dt; // Update angular acceration in body frame
    }
    /*Enforcing boundary conditions*/
    set_simple_bc(config, geometry, beam);

    if (check_energy_balance) {
        kinetic_energy_update(N, M, J_u, v_trans, v_rot, KE);
    }
}

void calc_initial_accelerations(const Config &config, const Geometry &geometry, const Borehole &borehole,
                                BeamSystem &beam) {
    calc_forces(config, geometry, borehole, beam);
    for (Index i = 0; i < geometry.get_N(); i++) {
        const Scalar M = beam.M[i];
        const Vec3 f = beam.R_ext_trans[i] - beam.R_int_trans[i];
        beam.a_trans[i] = f / M;

        const Mat3 U = beam.d_rot[i].to_matrix();
        const Mat3 &J_u = beam.J_u[i].asDiagonal();
        const Vec3 &m = beam.R_ext_rot[i] - beam.R_int_rot[i];
        const Vec3 &omega_u = beam.v_rot[i];
        beam.a_rot[i] = J_u.inverse() * (U.transpose() * m - omega_u.cross(J_u * omega_u));
    }
}

void calc_forces(const Config &config, const Geometry &geometry, const Borehole &borehole, BeamSystem &beam) {
    const Index N = geometry.get_N();

    zero_internal_and_set_static_forces(N, beam.R_static_trans, beam.R_static_rot, beam.R_int_trans, beam.R_int_rot,
                                        beam.R_ext_trans, beam.R_ext_rot);

    add_user_defined_external_forces(config, geometry, borehole, beam.d_trans, beam.d_rot, beam.R_ext_trans,
                                     beam.R_ext_rot);

    calc_inner_forces<0>(config, geometry, beam);
    calc_inner_forces<1>(config, geometry, beam);

    if (config.contact_enabled) {
        calc_hole_contact_forces(config, N, borehole.get_N_hole_elements(), borehole.get_x(), beam.hole_index,
                                 borehole.get_r_hole_element(), geometry.get_ro(), geometry.get_X(), beam.d_trans,
                                 beam.d_rot, beam.v_trans, beam.v_rot, beam.R_ext_trans, beam.R_ext_rot);
    }

    if (config.rayleigh_damping_enabled) {
        add_mass_proportional_rayleigh_damping(N, config.alpha_rayleigh, beam.M, beam.v_trans, beam.R_int_trans,
                                               beam.J_u, beam.d_rot, beam.v_rot, beam.R_int_rot);
    }
}

// /*Bending. Euler Bernoulli is used, symmetrical cross section
// The original matrix reads
// EI/L³ *[[ 12  6L  -12  6L ]
//         [ 6L  4L² -6L  2L²]
//         [-12 -6L   12 -6L ]
//         [ 6L  2L² -6L  4L²]]
// multiplied by [w1 theta1, w2 theta2],
// however, it can be simplified, since w1=w2=0,
// rows 1 and 3 can be skipped resulting in
// EI/L²**[[ 6   6 ]*[[theta1 ]
//         [ 4L  2L]  [theta2 ]]
//         [-6  -6 ]
//         [ 2L  4L]]
// */

// // const Scalar Fy1 = E * I / (l0 * l0) * (6 * theta_2A + 6 * theta_z2);
// // const Scalar Mz1 = E * I / (l0 * l0) * (4 * l0 * theta_z1 + 2 * l0 * theta_z2);
// // const Scalar Fy2 = -Fy1;
// // const Scalar Mz2 = E * I / (l0 * l0) * (2 * l0 * theta_z1 + 4 * l0 * theta_z2);

// const Scalar F2 = E * I / (l0 * l0) * (6 * theta_2l + 6 * theta_5l);
// const Scalar M2 = E * I / (l0 * l0) * (4 * l0 * theta_2l + 2 * l0 * theta_5l);
// const Scalar F5 = -F2;
// const Scalar M5 = E * I / (l0 * l0) * (2 * l0 * theta_2l + 4 * l0 * theta_5l);

// R_int_l[1] = F2;
// R_int_l[4] = M2;
// R_int_l[7] = F5;
// R_int_l[11] = M5;

// // const Scalar Fz1 = E * I / (l0 * l0) * (6 * theta_3A + 6 * theta_3B);
// // const Scalar My1 = E * I / (l0 * l0) * (4 * l0 * theta_3A + 2 * l0 * theta_3B);
// // const Scalar Fz2 = -Fz1;
// // const Scalar My2 = E * I / (l0 * l0) * (2 * l0 * theta_3A + 4 * l0 * theta_3B);

// const Scalar F3 = E * I / (l0 * l0) * (6 * theta_3l + 6 * theta_6l);
// const Scalar M3 = E * I / (l0 * l0) * (4 * l0 * theta_3l + 2 * l0 * theta_6l);
// const Scalar F6 = -F3;
// const Scalar M6 = E * I / (l0 * l0) * (2 * l0 * theta_3l + 4 * l0 * theta_6l);

// R_int_l[2] = F3;
// R_int_l[5] = M3;
// R_int_l[8] = F6;
// R_int_l[10] = M6;

// return R_int_l;
//}
//
static void zero_internal_and_set_static_forces(const Index N, const vector<Vec3> &R_static_trans,
                                                const vector<Vec3> &R_static_rot, vector<Vec3> &R_int_trans,
                                                vector<Vec3> &R_int_rot, vector<Vec3> &R_ext_trans,
                                                vector<Vec3> &R_ext_rot) {
#pragma omp parallel for
    for (Index i = 0; i < N; i++) {
        R_int_trans[i] = {0, 0, 0};
    }
#pragma omp parallel for
    for (Index i = 0; i < N; i++) {
        R_int_rot[i] = {0, 0, 0};
    }
#pragma omp parallel for
    for (Index i = 0; i < N; i++) {
        R_ext_trans[i] = R_static_trans[i];
    }
#pragma omp parallel for
    for (Index i = 0; i < N; i++) {
        R_ext_rot[i] = R_static_rot[i];
    }
}

template <Index i_first>
static void calc_inner_forces(const Config &config, const Geometry &geometry, BeamSystem &beam) {
    assert(i_first == 0 || i_first == 1);
    const Index N = geometry.get_N();
    const Index Ne = N - 1;
    const vector<Vec3> &X = geometry.get_X();
    const Scalar E = config.E;
    const Scalar G = config.get_G();
    const Scalar beta_rayleigh = config.beta_rayleigh;
    const CorotationalFormulation corotational_formulation = config.corotational_formulation;

#pragma omp parallel for
    for (Index ie = i_first; ie < Ne; ie += 2) {
        Vec12 R_int_e;
        Scalar A{}, I_2{}, I_3{}, J{};

        geometry.calc_cross_section_properties(ie, A, I_2, I_3, J);

        // I_2 = M_PI * powi<4>(0.1) / 64;
        // I_3 = I_2;

        // Scalar h_2 = 0.1;
        // Scalar h_3 = 0.05;

        // I_2 = h_2 * powi<3>(h_3) / 12;
        // I_3 = h_3 * powi<3>(h_2) / 12;

        switch (corotational_formulation) {
        default: {
            assert(corotational_formulation == CorotationalFormulation::CRISFIELD);
            R_int_e = CrisfieldBeam::calc_element_inner_forces(ie, X, beam.d_trans, beam.d_rot, E, G, I_2, I_3, A, J,
                                                               beta_rayleigh, beam.v_trans, beam.v_rot);
            break;
        }
        case CorotationalFormulation::BATTINI: {
            // have to flip these to get consistent behaviour compared to crisfield.. not exactly sure why
            const Scalar Iy = I_3;
            const Scalar Iz = I_2;
            R_int_e = BattiniBeam::global_internal_forces(ie, X, beam.d_trans, beam.d_rot, A, Iy, Iz, J, E, G);
            break;
        }
        }
        beam.R_int_trans[ie] += R_int_e.segment(0, 3);
        beam.R_int_rot[ie] += R_int_e.segment(3, 3);
        beam.R_int_trans[ie + 1] += R_int_e.segment(6, 3);
        beam.R_int_rot[ie + 1] += R_int_e.segment(9, 3);
    }

    /*------------------------------------------  # - R: [0,0,0, 0, 600000, 0]
    #   rel_loc: 1--------------------------
    To avoid race conditions the following pattern is used to compute and
    assemble the internal forces:

    Example shows a problem of 2 threads and 12 elements.

    first loop:
     e0,  e1,  e2,  e3,  e4,  e5,  e6,  e7,  e8,  e9, e10, e11
    [r0,  r0,  r0,  --,  --,  --,  r1,  r1,  r1,  --,  --, -- ]

    second loop
     e0,  e1,  e2,  e3,  e4,  e5,  e6,  e7,  e8,  e9, e10, e11
    [--,  --,  --,  r0,  r0,  r0,  --,  --,  --,  r1,  r1,  r1 ]

    ensuring that the memory region worked on by each thread is compact
    and, but that no two threads can update the same nodal force at the
    same time.
    --------------------------------------------------------------------*/
    // #pragma omp parallel
    //     {
    //         Index num_threads = omp_get_num_threads();
    //         Index bin_size = ceil((Scalar)Ne / (2 * num_threads));
    //         // cout << "bin sz " << bin_size << endl;
    //         assert(2 * bin_size * num_threads >= Ne);
    //         Index r = omp_get_thread_num();
    //         Index a = 2 * r * bin_size;
    //         assert(a < Ne);
    //         Index b = min(a + bin_size, Ne);
    //         // printf("first a %i,b %i\n", a, b);

    //         for (Index ie = a; ie < b; ie++)
    //         {
    //             calc_element_inner_forces(ie, X.data(), beam.d_trans.data(), beam.d_rot.data(),
    //                                       beam.R_int_trans.data(), beam.R_int_rot.data(),
    //                                       geometry.ri_e(ie), geometry.ro_e(ie), E, G);
    //         }
    //     }
    // #pragma omp parallel
    //     {
    //         Index num_threads = omp_get_num_threads();
    //         Index bin_size = ceil((Scalar)Ne / (2 * num_threads));
    //         assert(2 * bin_size * num_threads >= Ne);
    //         Index r = omp_get_thread_num();
    //         Index a = (2 * r + 1) * bin_size;
    //         assert(a < Ne);
    //         Index b = min(a + bin_size, Ne);
    //         // printf("second a %i,b %i\n", a, b);

    //         for (Index ie = a; ie < b; ie++)
    //         {
    //             calc_element_inner_forces(ie, X.data(), beam.d_trans.data(), beam.d_rot.data(),
    //                                       beam.R_int_trans.data(), beam.R_int_rot.data(),
    //                                       geometry.ri_e(ie), geometry.ro_e(ie), E, G);
    //         }
    //     }
}

static void velocity_update_partial_OLD(Scalar dt, Index N, const Scalar *__restrict__ M, const Vec3 *__restrict__ J_u,
                                        const Vec3 *__restrict__ R_int_trans, const Vec3 *__restrict__ R_int_rot,
                                        const Vec3 *__restrict__ R_ext_trans, const Vec3 *__restrict__ R_ext_rot,
                                        Vec3 *__restrict__ v_trans, Vec3 *__restrict__ v_rot) {
#pragma omp parallel
    {
#pragma omp for nowait
        for (Index i = 0; i < N; i++) {
            v_trans[i] += 0.5 * dt * (R_ext_trans[i] - R_int_trans[i]) / M[i];

            DEBUG_ONLY(if (i == 1) cout << "v_trans:\n" << v_trans[i] << endl;);
        }
#pragma omp for
        for (Index i = 0; i < N; i++) {
            /*omega_u_dot = J_u^-1 * (R_rot_u - S(omega_u)*J_u*omega_u)*/
            const Vec3 R_rot_u = R_ext_rot[i] - R_int_rot[i];

            DEBUG_ONLY(if (i == 1) cout << "R_ext_rot\n"
                                        << R_ext_rot[i] << endl
                                        << "R_int_rot\n"
                                        << R_int_rot[i] << endl
                                        << "R_rot\n"
                                        << R_rot_u << endl);
            Vec3 &omega_u = v_rot[i];

            Mat3 JJu = J_u[i].asDiagonal();
            Vec3 Ju = J_u[i];
            Vec3 rot_term_orig = omega_u.cross(Vec3{J_u[i].array() * omega_u.array()});

            Vec3 rot_term_new = skew(omega_u) * JJu * omega_u;

            cout << "rot_term_orig " << rot_term_orig << endl;
            cout << "rot_term_new " << rot_term_new << endl;

            Vec3 omega_u_dot_new;
            omega_u_dot_new.x() = (R_rot_u.x() - (Ju.z() - Ju.y()) * omega_u.y() * omega_u.z()) / Ju.x();
            omega_u_dot_new.y() = (R_rot_u.y() - (Ju.x() - Ju.z()) * omega_u.x() * omega_u.z()) / Ju.y();
            omega_u_dot_new.z() = (R_rot_u.z() - (Ju.y() - Ju.x()) * omega_u.x() * omega_u.y()) / Ju.z();

            const Vec3 omega_u_dot = (R_rot_u - rot_term_orig).array() / J_u[i].array();

            cout << "omega_u_dot\n" << omega_u_dot << endl;
            cout << "omega_u_dot_new\n" << omega_u_dot_new << endl;
            // cout << "diff: \n"
            //      << omega_u_dot - omega_u_dot_new << endl;

            // const Vec3 omega_u_dot = (R_rot_u - omega_u.cross(Vec3{J_u[i].array() * omega_u.array()})).array() /
            // J_u[i].array();
            //  Scalar tol = 0.00001;
            //  if (i != 0 && abs(omega_u[0]) > tol)
            //  {
            //      cout << "omega_u \n"
            //           << omega_u << endl;
            //      assert(false);
            //  }
            omega_u += 0.5 * dt * omega_u_dot;
            DEBUG_ONLY(if (i == 1) cout << "omega_u:\n" << v_rot[i] << endl;);
            // if (n_glob == 10000 && i > 10)
            //     omega_u[0] = 0.1;
        }
    }
}

static void displacement_update(Scalar dt, Index N, vector<Vec3> &v_trans, vector<Vec3> &v_rot, vector<Vec3> &d_trans,
                                vector<Quaternion> &d_rot) {
#pragma omp parallel for

    for (Index i = 0; i < N; i++) {
        d_trans[i] += dt * v_trans[i];
        DEBUG_ONLY(if (i == 1) cout << "d_trans:\n" << d_trans[i] << endl;);
    }
#pragma omp parallel for
    for (Index i = 0; i < N; i++) {
        // Mat3 U = d_rot[i].to_matrix();
        // U = U * (Mat3::Identity() - 0.5 * dt * skew(v_rot[i])).inverse() *
        //     (Mat3::Identity() + 0.5 * dt * skew(v_rot[i]));
        // assert(U.allFinite());
        // assert(is_orthogonal(U));
        // d_rot[i].from_matrix(U);
        // assert(is_close(d_rot[i].norm(), 1.0));

        Quaternion &q = d_rot[i];
        const Vec3 &omega_u = v_rot[i];

        const Vec3 omega = q.rotate_vector(
            omega_u); // perhaps possible to simplify this and not having to first convert omega to the global frame

        // const Vec3 omega = q.rotate_vector_reversed(omega_u); // perhaps possible to simplify this and not having
        // to first convert omega to the global frame
        if (i == 1) {
            cout << "i " << i << endl;
            cout << "omega " << omega << endl;
            cout << "omega mag " << omega.norm() << endl;
        }
        // assert(is_close(omega_u.y(), 0, 1e-4) && is_close(omega_u.z(), 0, 1e-4));
        assert(is_close(q.norm(), 1.0));
    }
}

static void calc_delta_d(Scalar dt, Index N, vector<Vec3> &delta_d_trans, vector<Vec3> &delta_d_rot,
                         const vector<Vec3> &v_trans, const vector<Vec3> &v_rot) {
#pragma omp parallel
    {
#pragma omp for nowait
        for (Index i = 0; i < N; i++) {
            delta_d_trans[i] = dt * v_trans[i];
        }
#pragma omp for
        for (Index i = 0; i < N; i++) {
            delta_d_rot[i] = dt * v_rot[i];
        }
    }
}
static void work_update_partial(Index N, const vector<Vec3> &delta_d_trans, const vector<Vec3> &delta_d_rot,
                                const vector<Vec3> &R_int_trans, const vector<Vec3> &R_int_rot,
                                const vector<Vec3> &R_ext_trans, const vector<Vec3> &R_ext_rot, Scalar &W_ext,
                                Scalar &W_int) {
#pragma omp parallel for reduction(+ : W_int) reduction(+ : W_ext)
    for (Index i = 0; i < N; i++) {
        /*The rotational dofs should have been rotated to the body frame allready*/
        W_int += 0.5 * (delta_d_trans[i].dot(R_int_trans[i]) + delta_d_rot[i].dot(R_int_rot[i]));
        W_ext += 0.5 * (delta_d_trans[i].dot(R_ext_trans[i]) + delta_d_rot[i].dot(R_ext_rot[i]));
    }
}

static void kinetic_energy_update(Index N, const vector<Scalar> &M, const vector<Vec3> &J_u,
                                  const vector<Vec3> &v_trans, const vector<Vec3> &v_rot, Scalar &KE) {
    KE = 0;
#pragma omp parallel for reduction(+ : KE)
    for (Index i = 0; i < N; i++) {
        const Vec3 &vt = v_trans[i];
        const Vec3 &omega_u = v_rot[i];
        KE += 0.5 * M[i] * vt.squaredNorm() + 0.5 * omega_u.dot(Vec3{J_u[i].array() * omega_u.array()});
    }
}

static void rotate_moment_to_body_frame(Index N, const vector<Quaternion> &d_rot, vector<Vec3> &R_int_rot,
                                        vector<Vec3> &R_ext_rot) {
#pragma omp parallel for
    for (Index i = 0; i < N; i++) {
        if (i == 1) {
            cout << "R_ext_rot before\n" << R_ext_rot[i] << endl;
        }
        //     Mat3 U = d_rot[i].to_matrix();
        //     R_int[i].rot = U.transpose() * R_int[i].rot;
        //     R_ext[i].rot = U.transpose() * R_ext[i].rot;
        const Quaternion &q = d_rot[i];

        R_int_rot[i] = q.rotate_vector_reversed(R_int_rot[i]);
        R_ext_rot[i] = q.rotate_vector_reversed(R_ext_rot[i]);

        if (i == 1) {
            cout << "R_ext_rot after\n" << R_ext_rot[i] << endl;
        }
    }
}

static void add_mass_proportional_rayleigh_damping(Index N, Scalar alpha, const vector<Scalar> &M,
                                                   const vector<Vec3> &v_trans, vector<Vec3> &R_int_trans,
                                                   const vector<Vec3> &J_u, const vector<Quaternion> &d_rot,
                                                   const vector<Vec3> &v_rot, vector<Vec3> &R_int_rot) {
#pragma omp parallel for
    for (Index i = 0; i < N; i++) {
        R_int_trans[i] += alpha * M[i] * v_trans[i];
    }
    // Test with rotational dofs also?
    // #pragma omp parallel for
    //     for (Index i = 0; i < N; i++)
    //     {
    //         Scalar alpha_rot = 0.3 * alpha;
    //         Vec3 R_damp_rot = alpha_rot * J_u[i].array() * v_rot[i].array();
    //         R_damp_rot = d_rot->rotate_vector(R_damp_rot);
    //         // if (i == 10 && n_glob % 100 == 0)
    //         // {
    //         //     cout << "M_w_1 prior: " << R_int_rot[i].x() << endl;
    //         //     cout << "M damp w_1: " << R_damp_rot.x() << endl;
    //         // }
    //         // cout << "R_int_rot \n"
    //         //      << R_int_rot[i] << endl
    //         //      << "R_damp\n"
    //         //      << R_damp_rot << endl;
    //         R_int_rot[i] += R_damp_rot;
    //     }
}

static void set_simple_bc(const Config &config, const Geometry &geometry, BeamSystem &beam) {
    Index N = geometry.get_N();

    vector<Vec3> &d_trans = beam.d_trans;
    vector<Quaternion> &d_rot = beam.d_rot;
    vector<Vec3> &v_trans = beam.v_trans;
    vector<Vec3> &v_rot = beam.v_rot;
    vector<Vec3> &delta_d_trans = beam.delta_d_trans; // design withoug having to set these later
    vector<Vec3> &delta_d_rot = beam.delta_d_rot;

    assert(N == d_trans.size() && N == d_rot.size() && N == v_trans.size() && N == v_rot.size());

    switch (config.bc_case) {
    case BC_Case::NONE:
        return;
    case BC_Case::CANTILEVER:
        d_trans[0] = {0, 0, 0};
        v_trans[0] = {0, 0, 0};
        delta_d_trans[0] = {0, 0, 0};
        // d_rot[0].from_matrix();
        d_rot[0] = config.bc_orientation_base;
        v_rot[0] = {0, 0, 0};
        delta_d_rot[0] = {0, 0, 0};
        break;
    case BC_Case::SIMPLY_SUPPORTED:
        d_trans[0] = {0, 0, 0};
        v_trans[0] = {0, 0, 0};
        d_trans[N - 1].y() = 0;
        d_trans[N - 1].z() = 0;
        v_trans[N - 1].y() = 0;
        v_trans[N - 1].z() = 0;
        assert(false); // fix delta d trans and rot
        break;
    default:
        assert(false);
    }
}