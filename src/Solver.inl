
#pragma once
#include "../include/Solver.hpp"



inline void calc_element_inner_forces(Index ie, const Vec3 *__restrict__ X, const Vec3 *__restrict__ d_trans,
                                      const Quaternion *__restrict__ d_rot, Vec3 *__restrict__ R_int_trans,
                                      Vec3 *__restrict__ R_int_rot, Scalar ri_e, Scalar ro_e, Scalar E, Scalar G)
{

    // if (n_glob == 2)
    // {
    //     int a = 1;
    // }

    // assert(i&e < X.size() - 1);
    const Vec3 &X1 = X[ie];
    const Vec3 &X2 = X[ie + 1];
    const Vec3 &d1 = d_trans[ie];
    const Vec3 &d2 = d_trans[ie + 1];

    // cout << "d1 " << d1.transpose() << endl;
    // cout << "d2 " << d2.transpose() << endl;

    const Mat3 T = d_rot[ie].to_matrix();
    const Vec3 &t1 = T.col(0);
    const Vec3 &t2 = T.col(1);
    const Vec3 &t3 = T.col(2);

    const Mat3 U = d_rot[ie + 1].to_matrix();
    const Vec3 &u1 = U.col(0);
    const Vec3 &u2 = U.col(1);
    const Vec3 &u3 = U.col(2);

    // cout << "U from quat\n"
    //      << U << endl;

    // cout << "u1 u2 u3\n"
    //      << u1.transpose() << ",\n"
    //      << u2.transpose() << ",\n"
    //      << u3.transpose() << endl;

    // cout << "t1 t2 t3\n"
    //      << t1.transpose() << ",\n"
    //      << t2.transpose() << ",\n"
    //      << t3.transpose() << endl;

    // calculate the first unit vector
    const Scalar l0 = (X2 - X1).norm();
    const Scalar ln = (X2 + d2 - (X1 + d1)).norm();

    const Vec3 e1 = (X2 + d2 - (X1 + d1)) / ln;
    assert(is_close(e1.norm(), 1.0));

    // cout << "e1 = " << e1.transpose() << endl;
    // Calculate intermediate rotation between triads
    const Mat3 DeltaR = U * T.transpose();
    Quaternion qDeltaR;
    qDeltaR.from_matrix(DeltaR);
    const Vec3 gamma_half = qDeltaR.q / qDeltaR.q0; // eq 16.34 divided by 2
    const Mat3 S = skew_symmetric(gamma_half);
    const Mat3 DeltaR_m = Mat3::Identity() + 1 / (1 + 0.25 * gamma_half.dot(gamma_half)) * (S + 0.5 * S * S);
    const Mat3 R_ = DeltaR_m * T;
    const Vec3 &r1 = R_.col(0);
    const Vec3 &r2 = R_.col(1);
    const Vec3 &r3 = R_.col(2);
    // cout << "R_\n"
    //      << R_ << endl;

    // computing element unit vectors e2 and e3 by rotating the unit vector
    // r1 on to e1.
    const Vec3 e2 = r2 - r2.dot(e1) / (1 + r1.dot(e1)) * (e1 + r1);
    const Vec3 e3 = r3 - r3.dot(e1) / (1 + r1.dot(e1)) * (e1 + r1);
    // const Vec3 e3 = e1.cross(e2);
    // cout << "e2 = " << e2.transpose() << endl;
    // cout << "e3 = " << e3.transpose() << endl;

    // compute local displacements from global displacements.
    // Taking u_l = ln-l0 is not recommended since subtracting two large
    // numbers may be inaccurate with limited precision. Better to adopt
    // ul = ln - l0 = (ln - l0)*(ln + l0)/(ln + l0) = (ln^2 + l0^2)/(ln + l0)

    const Scalar ul = (ln * ln - l0 * l0) / (ln + l0);
    assert(is_close(ul, ln - l0));

    const Scalar theta_l1 = asin(0.5 * (-t3.dot(e2) + t2.dot(e3)));
    const Scalar theta_l2 = asin(0.5 * (-t2.dot(e1) + t1.dot(e2)));
    const Scalar theta_l3 = asin(0.5 * (-t3.dot(e1) + t1.dot(e3)));
    const Scalar theta_l4 = asin(0.5 * (-u3.dot(e2) + u2.dot(e3)));
    const Scalar theta_l5 = asin(0.5 * (-u2.dot(e1) + u1.dot(e2)));
    const Scalar theta_l6 = asin(0.5 * (-u3.dot(e1) + u1.dot(e3)));
    // cout << "n_glob " << n_glob << endl;
    // assert(theta_l1 == 0);
    // assert(theta_l4 == 0);

#define MAX_ANGLE 5 * M_PI / 180
    assert(abs(theta_l1) < MAX_ANGLE);
    assert(abs(theta_l2) < MAX_ANGLE);
    assert(abs(theta_l3) < MAX_ANGLE);
    assert(abs(theta_l4) < MAX_ANGLE);
    assert(abs(theta_l5) < MAX_ANGLE);
    assert(abs(theta_l6) < MAX_ANGLE);
#undef MAX_ANGLE

    Eigen::Matrix<Scalar, 12, 7> F_transpose; // Transpose of F where zero rows in F are excluded
    constexpr Index f4 = 0;
    constexpr Index f5 = 1;
    constexpr Index f6 = 2;
    constexpr Index f7 = 3;
    constexpr Index f10 = 4;
    constexpr Index f11 = 5;
    constexpr Index f12 = 6;

    // calculate F connecting infinitesimal global and local variable
    // optimize zero rows later

    // Vec12 f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12;
    // f1.setConstant(0);
    // f2.setConstant(0);
    // f3.setConstant(0);
    // f8.setConstant(0);
    // f9.setConstant(0);

    F_transpose.col(f7) << -e1, Vec3::Zero(), e1, Vec3::Zero(); //(17.19)

    const Mat3 A = 1.0 / ln * (Mat3::Identity() - e1 * e1.transpose());
    assert(is_close((A - A.transpose()).norm(), 0.0));

    Mat3 L1r2 = 0.5 * r2.dot(e1) * A + 0.5 * A * r2 * (e1 + r1).transpose();
    Mat3 L1r3 = 0.5 * r3.dot(e1) * A + 0.5 * A * r3 * (e1 + r1).transpose();
    Mat3 L2r2 = 0.5 * skew_symmetric(r2) - 0.25 * r2.transpose() * e1 * skew_symmetric(r1) -
                0.25 * skew_symmetric(r2) * e1 * (e1 + r1).transpose();
    Mat3 L2r3 = 0.5 * skew_symmetric(r3) - 0.25 * r3.transpose() * e1 * skew_symmetric(r1) -
                0.25 * skew_symmetric(r3) * e1 * (e1 + r1).transpose();

    Eigen::Matrix<Scalar, 12, 3> Lr2, Lr3;
    Lr2 << L1r2, L2r2, -L1r2, L2r2;
    Lr3 << L1r3, L2r3, -L1r3, L2r3;

    Vec12 h1, h2, h3, h4, h5, h6;
    h1 << Vec3::Zero(), -t3.cross(e2) + t2.cross(e3), Vec3::Zero(), Vec3::Zero();
    h2 << A * t2, -t2.cross(e1) + t1.cross(e2), -A * t2, Vec3::Zero();
    h3 << A * t3, -t3.cross(e1) + t1.cross(e3), -A * t3, Vec3::Zero();
    h4 << Vec3::Zero(), Vec3::Zero(), Vec3::Zero(), -u3.cross(e2) + u2.cross(e3);
    h5 << A * u2, Vec3::Zero(), -A * u2, -u2.cross(e1) + u1.cross(e2);
    h6 << A * u3, Vec3::Zero(), -A * u3, -u3.cross(e1) + u1.cross(e3);

    // cout << "h3 " << h3.transpose() << endl;
    // cout << "h6 " << h6.transpose() << endl;

    F_transpose.col(f4) = 1 / (2 * cos(theta_l1)) * (Lr3 * t2 - Lr2 * t3 + h1);
    F_transpose.col(f5) = 1 / (2 * cos(theta_l2)) * (Lr2 * t1 + h2);
    F_transpose.col(f6) = 1 / (2 * cos(theta_l3)) * (Lr3 * t1 + h3);
    F_transpose.col(f10) = 1 / (2 * cos(theta_l4)) * (Lr3 * u2 - Lr2 * u3 + h4);
    F_transpose.col(f11) = 1 / (2 * cos(theta_l5)) * (Lr2 * u1 + h5); // seems to be an error in the book for f11 and f12. it should be + in front of h6 and h6 (it's + in the paper)
    F_transpose.col(f12) = 1 / (2 * cos(theta_l6)) * (Lr3 * u1 + h6);

    // cout << "f6\n"
    //      << f6 << "\nf12\n"
    //      << f12 << endl;

    // Eigen::Matrix<Scalar, 12, 12> F;
    // F << f1.transpose(),
    //     f2.transpose(),
    //     f3.transpose(),
    //     f4.transpose(),
    //     f5.transpose(),
    //     f6.transpose(),
    //     f7.transpose(),
    //     f8.transpose(),
    //     f9.transpose(),
    //     f10.transpose(),
    //     f11.transpose(),
    //     f12.transpose();

    // cout << "F^T\n"
    //      << F.transpose() << endl;

    /*Calculate local internal forces based on linear 3D beam theory*/
    const Vec7 R_int_e_l = calc_element_forces_local(ri_e, ro_e, l0, E, G, ul,
                                                     theta_l1, theta_l2, theta_l3, theta_l4, theta_l5, theta_l6);

    // cout << "R_int_l:\n"
    //      << R_int_l << endl;

    const Vec12 R_int_e = F_transpose * R_int_e_l;
    assert(R_int_e.allFinite());
    // cout << "R_int:\n"
    //      << R_int << endl;

    // Vec3 f1, m1;
    // Vec3 f2, m2;

    // calc_element_forces_local_TEST(ri_e, ro_e, l0, E, G, ul,
    //                                theta_l1, theta_l2, theta_l3, theta_l4, theta_l5, theta_l6, f1, m1, f2, m2);
    // Mat3 Etri;
    // Etri.col(0) = e1;
    // Etri.col(1) = e2;
    // Etri.col(2) = e3;
    // f1 = Etri * f1;
    // f2 = Etri * f2;
    // m1 = Etri * m1;
    // m2 = Etri * m2;

    // Vec12 R_rotated_eng;
    // R_rotated_eng << f1, m1, f2, m2;

    // cout << "R_int_e\n"
    //      << R_int_e << endl;
    // cout << "R_rotated_eng\n"
    //      << R_rotated_eng << endl;

    // Vec12 rel_diff_percent = 100 * (R_int_e - R_rotated_eng).array() / R_int_e.array();
    // for (Index i = 0; i < 12; i++)
    // {
    //     if (!isnan(rel_diff_percent[i]) && R_int_e[i] > 100)
    //     {
    //         // printf("rel diff too big\n");

    //         assert(abs(rel_diff_percent[i] < 1));
    //     }
    // }
    // R_int_trans[ie] += f1;
    // R_int_rot[ie] += m1;
    // R_int_trans[ie + 1] += f2;
    // R_int_rot[ie + 1] += m2;
    // cout << "rel diff %\n"
    //      << rel_diff_percent << endl;
    /*R is the residual force moved to the right hand side. If we have:
    M*a + R_int = R_ext, this is here instead handled as:
    M*a = R, meaning that R_int has negative contributions to R*/
    R_int_trans[ie] += R_int_e.segment(0, 3);
    R_int_rot[ie] += R_int_e.segment(3, 3);
    R_int_trans[ie + 1] += R_int_e.segment(6, 3);
    R_int_rot[ie + 1] += R_int_e.segment(9, 3);
}

inline Vec7 calc_element_forces_local(Scalar ri, Scalar ro, Scalar l0, Scalar E, Scalar G, Scalar ul,
                                      Scalar theta_1l, Scalar theta_2l, Scalar theta_3l, Scalar theta_4l,
                                      Scalar theta_5l, Scalar theta_6l)
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

    /*Bending: Euler bernoulli with only angle dofs:
    k = EI/L [[4 2],
              [2 4]]
    */
    const Scalar M2 = E * I / l0 * (4 * theta_2l + 2 * theta_5l);
    const Scalar M5 = E * I / l0 * (2 * theta_2l + 4 * theta_5l);

    const Scalar M3 = E * I / l0 * (4 * theta_3l + 2 * theta_6l);
    const Scalar M6 = E * I / l0 * (2 * theta_3l + 4 * theta_6l);

    // Vec12 R_int_l = {F1, 0, 0, M1, M2, M3, F4, 0, 0, M4, M5, M6};
    const Vec7 R_int_e_l = {M1, M2, M3, F4, M4, M5, M6};
    // cout << "fix\n";
    // //Vec12 R_int_l = {0, 0, 0, 0, M2, M3, 0, 0, 0, 0, M5, M6};
    // Vec12 R_int_l = {0, 0, 0, 0, M2, M3, 0, 0, 0, 0, M5, M6};
    assert(R_int_e_l.allFinite());
    return R_int_e_l;
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

inline void assemble(const Config &config, const Geometry &geometry, BeamSystem &beam_system)
{
    const Index N = geometry.get_N();
    const Index Ne = N - 1;

    const vector<Vec3> &X = geometry.get_X();
    const Scalar E = config.E;
    const Scalar G = config.get_G();

/*Start by setting internal forces to zero and external forces to the static loads*/
#pragma omp parallel for
    for (Index i = 0; i < N; i++)
    {
        beam_system.R_int_trans[i] = {0, 0, 0};
        beam_system.R_int_rot[i] = {0, 0, 0};
        beam_system.R_ext_trans[i] = beam_system.R_static_trans[i];
        beam_system.R_ext_rot[i] = beam_system.R_static_rot[i];
    }

#pragma omp parallel for
    for (Index ie = 0; ie < Ne; ie++)
    {
        if (ie % 2 == 0)
        {
            calc_element_inner_forces(ie, X.data(), beam_system.d_trans.data(), beam_system.d_rot.data(),
                                      beam_system.R_int_trans.data(), beam_system.R_int_rot.data(),
                                      geometry.ri_e(ie), geometry.ro_e(ie), E, G);
        }
    }
#pragma omp parallel for
    for (Index ie = 0; ie < Ne; ie++)
    {

        if (ie % 2 == 1)
        {
            calc_element_inner_forces(ie, X.data(), beam_system.d_trans.data(), beam_system.d_rot.data(),
                                      beam_system.R_int_trans.data(), beam_system.R_int_rot.data(),
                                      geometry.ri_e(ie), geometry.ro_e(ie), E, G);
        }
    }

    /*--------------------------------------------------------------------
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
    //             calc_element_inner_forces(ie, X.data(), beam_system.d_trans.data(), beam_system.d_rot.data(),
    //                                       beam_system.R_int_trans.data(), beam_system.R_int_rot.data(),
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
    //             calc_element_inner_forces(ie, X.data(), beam_system.d_trans.data(), beam_system.d_rot.data(),
    //                                       beam_system.R_int_trans.data(), beam_system.R_int_rot.data(),
    //                                       geometry.ri_e(ie), geometry.ro_e(ie), E, G);
    //         }
    //     }
}

inline void velocity_update_partial(Scalar dt, Index N, const Scalar *__restrict__ M, const Vec3 *__restrict__ J_u,
                                    const Vec3 *__restrict__ R_int_trans, const Vec3 *__restrict__ R_int_rot,
                                    const Vec3 *__restrict__ R_ext_trans, const Vec3 *__restrict__ R_ext_rot,
                                    Vec3 *__restrict__ v_trans, Vec3 *__restrict__ v_rot)
{
#pragma omp parallel
    {
#pragma omp for nowait
        for (Index i = 0; i < N; i++)
        {
            v_trans[i] += 0.5 * dt * (R_ext_trans[i] - R_int_trans[i]) / M[i];
        }
#pragma omp for
        for (Index i = 0; i < N; i++)
        {
            /*omega_u_dot = J_u^-1 * (R_rot_u - S(omega_u)*J_u*omega_u)*/
            const Vec3 R_rot_u = R_ext_rot[i] - R_int_rot[i];
            const Vec3 &omega_u = v_rot[i];
            const Vec3 omega_u_dot = (R_rot_u - omega_u.cross(Vec3{J_u[i].array() * omega_u.array()})).array() / J_u[i].array();
            v_rot[i] += 0.5 * dt * omega_u_dot;
        }
    }
}

inline void displacement_update(Scalar dt, Index N, Vec3 *__restrict__ v_trans, Vec3 *__restrict__ v_rot, Vec3 *__restrict__ d_trans,
                                Quaternion *__restrict__ d_rot)
{
#pragma omp parallel
    {
#pragma omp for nowait
        for (Index i = 0; i < N; i++)
        {
            d_trans[i] += dt * v_trans[i];
        }
#pragma omp for
        for (Index i = 0; i < N; i++)
        {
            // Mat3 U = d_rot[i].to_matrix();
            // U = U * (Mat3::Identity() - 0.5 * dt * skew_symmetric(v[i].rot)).inverse() *
            //     (Mat3::Identity() + 0.5 * dt * skew_symmetric(v[i].rot));
            // assert(U.allFinite());
            // d_rot[i].from_matrix(U);

            Quaternion &q = d_rot[i];
            const Vec3 &omega_u = v_rot[i];
            const Vec3 omega = q.rotate_vector(omega_u); // perhaps possible to simplify this and not having to first convert omega to the global frame
            q.compound_rotate(dt * omega);
        }
    }
}

inline void calc_delta_d(Scalar dt, Index N, Vec3 *__restrict__ delta_d_trans, Vec3 *__restrict__ delta_d_rot,
                         const Vec3 *__restrict__ v_trans, const Vec3 *__restrict__ v_rot)
{
#pragma omp parallel
    {
#pragma omp for nowait
        for (Index i = 0; i < N; i++)
        {
            delta_d_trans[i] = dt * v_trans[i];
        }
#pragma omp for
        for (Index i = 0; i < N; i++)
        {
            delta_d_rot[i] = dt * v_rot[i];
        }
    }
}
inline void work_update_partial(Index N, const Vec3 *__restrict__ delta_d_trans, const Vec3 *__restrict__ delta_d_rot,
                                const Vec3 *__restrict__ R_int_trans, const Vec3 *__restrict__ R_int_rot,
                                const Vec3 *__restrict__ R_ext_trans, const Vec3 *__restrict__ R_ext_rot,
                                Scalar &W_ext, Scalar &W_int)
{
#pragma omp parallel for reduction(+ : W_int) reduction(+ : W_ext)
    for (Index i = 0; i < N; i++)
    {
        /*The rotational dofs should have been rotated to the body frame allready*/
        W_int += 0.5 * (delta_d_trans[i].dot(R_int_trans[i]) + delta_d_rot[i].dot(R_int_rot[i]));
        W_ext += 0.5 * (delta_d_trans[i].dot(R_ext_trans[i]) + delta_d_rot[i].dot(R_ext_rot[i]));
    }
}

inline void kinetic_energy_update(Index N, const Scalar *__restrict__ M, const Vec3 *__restrict__ J_u,
                                  const Vec3 *__restrict__ v_trans, const Vec3 *__restrict__ v_rot, Scalar &KE)
{
    KE = 0;
#pragma omp parallel for reduction(+ : KE)
    for (Index i = 0; i < N; i++)
    {
        const Vec3 &vt = v_trans[i];
        const Vec3 &omega_u = v_rot[i];
        KE += 0.5 * M[i] * vt.squaredNorm() + 0.5 * omega_u.dot(Vec3{J_u[i].array() * omega_u.array()});
    }
}

inline void rotate_moment_to_body_frame(Index N, const Quaternion *__restrict__ d_rot,
                                        Vec3 *__restrict__ R_int_rot, Vec3 *__restrict__ R_ext_rot)
{
#pragma omp parallel for
    for (Index i = 0; i < N; i++)
    {
        //     Mat3 U = d_rot[i].to_matrix();
        //     R_int[i].rot = U.transpose() * R_int[i].rot;
        //     R_ext[i].rot = U.transpose() * R_ext[i].rot;
        const Quaternion &q = d_rot[i];
        R_int_rot[i] = q.rotate_vector_reversed(R_int_rot[i]);
        R_ext_rot[i] = q.rotate_vector_reversed(R_ext_rot[i]);
    }
}

inline void add_mass_proportional_rayleigh_damping(Index N, Scalar alpha, const Scalar *__restrict__ M,
                                                   const Vec3 *__restrict__ v_trans, Vec3 *__restrict__ R_int_trans)
{
#pragma omp parallel for
    for (Index i = 0; i < N; i++)
    {
        R_int_trans[i] += alpha * M[i] * v_trans[i];
    }
}

// inline void step_central_differences(Scalar dt, Index N, Vec3Quat *__restrict__ u, Vec3Vec3 *__restrict__ v,
//                                      const Scalar *__restrict__ M, const Vec3 *__restrict__ J_u,
//                                      const Vec3Vec3 *__restrict__ R_int, const Vec3Vec3 *__restrict__ R_ext,
//                                      bool check_energy_balance, Scalar &W_int, Scalar &W_ext, Scalar &KE)
// {

//     // print_std_vector(R, "R");

//     // print_std_vector(R_static, "R_static");

//     // Vec3Quat::print_array(u, "u", true, true);

//     // print_std_vector(v, "v");

//     Scalar dW_int = 0, dW_ext = 0;
//     KE = 0;

// #pragma omp parallel for reduction(+ : dW_int) reduction(+ : dW_ext) reduction(+ : KE)
//     for (Index i = 0; i < N; i++)
//     {
//         // cout << "remove!\n";
//         // if (i == 0)
//         //     continue;

//         const Vec3 &R_int_trans = R_int[i].trans;
//         const Vec3 &R_int_rot = R_int[i].rot;
//         const Vec3 &R_ext_trans = R_ext[i].trans;
//         const Vec3 &R_ext_rot = R_ext[i].rot;
//         assert(R_int_trans.allFinite());
//         assert(R_int_rot.allFinite());
//         assert(R_ext_trans.allFinite());
//         assert(R_ext_rot.allFinite());

//         /*--------------------------------------------------------------------
//         Translations
//         --------------------------------------------------------------------*/
//         v[i].trans += dt * (R_ext_trans - R_int_trans) / M[i];
//         u[i].trans += dt * v[i].trans;
//         // cout << "v " << v[i].trans << endl
//         //      << "u " << u[i].trans << endl;
//         /*--------------------------------------------------------------------
//         Rotations
//         --------------------------------------------------------------------*/
//         Quaternion &q = u[i].rot;
//         // Optimize this! use Quaternion product directly instead to rotate
//         Mat3 U = q.to_matrix();
//         // cout << "U prior\n"
//         //      << U << endl;
//         const Vec3 R_rot_u = U.transpose() * (R_ext_rot - R_int_rot); /*Transform nodal force to node frame*/
//         Vec3 &omega_u = v[i].rot;                                     /*The angular velocities are related to the body frame*/
//         const Vec3 &J_u_i = J_u[i];

//         /*omega_u_dot = J_u^-1 * (R_rot_u - S(omega_u)*J_u*omega_u)*/
//         // cout << "R_rot_u " << R_rot_u << endl;
//         // cout << "J_u " << J_u_i << endl;
//         const Vec3 omega_u_dot = (R_rot_u - omega_u.cross(Vec3{J_u_i.array() * omega_u.array()})).array() / J_u_i.array();
//         // cout << "omega_u_dot " << omega_u_dot << endl;
//         omega_u += dt * omega_u_dot;
//         // cout << "omega_u " << omega_u << endl;
//         /*How to update rotations: For now I will use Hughes-Winges.
//         I could also consider using the Quaternion diff eqn to do it, and normalize it,
//         would probably be better/faster

//         Hughes Winges: (Computing inverse or direct solver?) note that thus differs
//         from the equation in the hopperstad lecture notes since here omega_u is
//         used instead of omega*/
//         U = U * (Mat3::Identity() - 0.5 * dt * skew_symmetric(omega_u)).inverse() *
//             (Mat3::Identity() + 0.5 * dt * skew_symmetric(omega_u));
//         assert(U.allFinite());

//         /*Quaternion compound rotation*/
//         // const Vec3 omega = U * omega_u;
//         // const Vec3 Delta_Theta_pseudo = dt * omega;
//         // const Scalar Delta_Theta_val = Delta_Theta_pseudo.norm();
//         // Quaternion delta_q;
//         // delta_q.q0 = cos(Delta_Theta_pseudo);
//         // delta_q.q1 =
//         //     quat = quat;

//         // Eigen::Quaternion<Scalar> qu{U};
//         // Vec3 a;
//         // Vec3 b = qu*a;

//         // U = U * (Mat3::Identity() + dt * skew_symmetric(omega_u));
//         //   if (i == N - 1)
//         //   {
//         //       cout << "omega_u " << omega_u.transpose() << endl;
//         //   }
//         //   cout << "U\n"
//         //        << U << endl;
//         u[i].rot.from_matrix(U);

//         // R[i].set_zero(); /*Setting the load vector to zero so it's ready for assembly in the next timestep*/
//     }
// }

inline void set_simple_bc(BC_Case bc_case, const Geometry &geometry, BeamSystem &beam_system)
{
    Index N = geometry.get_N();

    vector<Vec3> d_trans = beam_system.d_trans;
    vector<Quaternion> &d_rot = beam_system.d_rot;
    vector<Vec3> &v_trans = beam_system.v_trans;
    vector<Vec3> &v_rot = beam_system.v_rot;

    assert(N == d_trans.size() && N == d_rot.size() && N == v_trans.size() && N == v_rot.size());

    switch (bc_case)
    {
    case BC_Case::NONE:
        return;
    case BC_Case::CANTILEVER:
        d_trans[0] = {0, 0, 0};
        v_trans[0] = {0, 0, 0};
        d_rot[0].from_matrix(Mat3::Identity());
        v_rot[0] = {0, 0, 0};
        break;
    case BC_Case::SIMPLY_SUPPORTED:
        d_trans[0] = {0, 0, 0};
        v_trans[0] = {0, 0, 0};
        d_trans[N - 1] = {0, 0, 0};
        v_trans[N - 1] = {0, 0, 0};
        break;
    default:
        assert(false);
    }
}