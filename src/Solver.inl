
#pragma once
#include "../include/Solver.hpp"

#define calc_element_inner_forces calc_element_inner_forces_crisfield

inline void calc_element_inner_forces_battini(Index ie, const Vec3 *__restrict__ X, const Vec3 *__restrict__ d_trans,
                                              const Quaternion *__restrict__ d_rot, Vec3 *__restrict__ R_int_trans,
                                              Vec3 *__restrict__ R_int_rot, Scalar ri_e, Scalar ro_e, Scalar youngs, Scalar G)
{
    DEBUG_ONLY(cout << "!!!!!!!!!!!!BATTINI!!!!!!!!!!\n");
    const Vec3 &X1 = X[ie];
    const Vec3 &X2 = X[ie + 1];
    const Vec3 &u1 = d_trans[ie];
    const Vec3 &u2 = d_trans[ie + 1];
    const Mat3 Rg1 = d_rot[ie].to_matrix();
    const Mat3 Rg2 = d_rot[ie + 1].to_matrix();

    const Scalar ln = (X2 + u2 - X1 - u1).norm();
    const Scalar l0 = (X2 - X1).norm();

    /*Create element frame*/
    const Vec3 e1 = (X2 + u2 - X1 - u1) / ln;
    const Vec3 q1 = Rg1 * Vec3{0, 1, 0};
    const Vec3 q2 = Rg2 * Vec3{0, 1, 0};
    const Vec3 q = 0.5 * (q1 + q2);
    const Vec3 e3 = (e1.cross(q)).normalized();
    const Vec3 e2 = e3.cross(e1);
    assert(is_close(e1.norm(), 1.0));
    assert(is_close(e3.norm(), 1.0) && is_close(e3.dot(e1), 0.0));
    assert(is_close(e2.norm(), 1.0) && is_close(e2.dot(e1), 0.0));
    Mat3 Rr;
    Rr << e1, e2, e3;
    DEBUG_ONLY(
        cout << "e1 = " << e1.transpose() << endl;
        cout << "e2 = " << e2.transpose() << endl;
        cout << "e3 = " << e3.transpose() << endl;);

    const Scalar ul = ln - l0;

    const Mat3 Rl1 = Rr.transpose() * Rg1;
    const Mat3 Rl2 = Rr.transpose() * Rg2;

    // Goal is to compute log(R).
    // we have R = exp(S(theta)), so log(R) = S(theta)?
    // I think I should use spurriers algorithm to extract a quaternion from R. And then extract the pseudo vector?
    Quaternion ql1, ql2;
    ql1.from_matrix(Rl1);
    ql2.from_matrix(Rl2);
    // omega = 2*q/q0
    const Vec3 theta_l1 = 2 * ql1.q / ql1.q0; // remember that omega = 2*tan(theta/2) is a very good approximation to theta for small angles
    const Vec3 theta_l2 = 2 * ql2.q / ql2.q0; // remember that omega = 2*tan(theta/2) is a very good approximation to theta for small angles
    constexpr Scalar THETA_MAX = 10 * M_PI / 180;
    assert(theta_l1.norm() < THETA_MAX);
    assert(theta_l2.norm() < THETA_MAX);
    DEBUG_ONLY(
        cerr << "theta_l1 = " << theta_l1.transpose() << endl;
        cerr << "theta_l2 = " << theta_l2.transpose() << endl;);

#ifndef NDEBUG
    // checking if the above corresponds with alternate definition (16.8 crisfield) allowed for 0 < |theta| <pi
    // checking for the first node:
    {
        Mat3 R = Rl1;
        Vec3 sin_theta_e = 0.5 * Vec3(R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1));
        Vec3 theta = sin_theta_e; // sin(theta)/theta are close for small angles
        cerr << "theta check " << theta << endl;
        assert(theta.isApprox(theta_l1, 1e-2));
    }
#endif

    // inverse transformation matrix
    Mat3 Ts_inv1, Ts_inv2;
    {
        constexpr Scalar tol = 1e-4;
        Vec3 theta_pseudo = theta_l1;
        Scalar theta = theta_pseudo.norm();
        Vec3 e = theta_pseudo.normalized();
        if (abs(theta) > tol)
        {
            Ts_inv1 = (theta / 2) / tan(theta / 2) * Mat3::Identity() + (1 - (theta / 2) / tan(theta / 2)) * e * e.transpose() - 0.5 * skew_symmetric(theta_pseudo);
        }
        else
        {
            Ts_inv1 = Mat3::Identity();
        }
        theta_pseudo = theta_l1;
        theta = theta_pseudo.norm();
        e = theta_pseudo.normalized();
        if (abs(theta) > tol)
        {
            Ts_inv2 = (theta / 2) / tan(theta / 2) * Mat3::Identity() + (1 - (theta / 2) / tan(theta / 2)) * e * e.transpose() - 0.5 * skew_symmetric(theta_pseudo);
        }
        else
        {
            Ts_inv2 = Mat3::Identity();
        }
    }
    DEBUG_ONLY(cout << "Ts_inv1\n"
                    << Ts_inv1 << endl
                    << "Ts_inv2\n"
                    << Ts_inv2 << endl;)

    using Mat7 = Eigen::Matrix<Scalar, 7, 7>;
    Mat7 Bl{Mat7::Zero()};

    Bl(0, 0) = 1;
    Bl.block<3, 3>(1, 1) = Ts_inv1;
    Bl.block<3, 3>(4, 4) = Ts_inv2;
    DEBUG_ONLY(cout << "Bl\n"
                    << Bl << endl;)

    using Vec12 = Eigen::Vector<Scalar, 12>;
    Vec12 r_T{Vec12::Zero()};
    r_T << -e1, Vec3::Zero(), e1, Vec3::Zero();
    using Mat12 = Eigen::Matrix<Scalar, 12, 12>;
    Mat12 E{Mat12::Zero()};
    E.block<3, 3>(0, 0) = Rr;
    E.block<3, 3>(3, 3) = Rr;
    E.block<3, 3>(6, 6) = Rr;
    E.block<3, 3>(9, 9) = Rr;

    Eigen::Matrix<Scalar, 3, 12> G_T;
    G_T.setZero();
    {
        const Vec3 _q = Rr.transpose() * q;
        const Vec3 _q1 = Rr.transpose() * q1;
        const Vec3 _q2 = Rr.transpose() * q2;

        const Scalar eta = _q.x() / _q.y();
        const Scalar eta11 = _q1.x() / _q.y();
        const Scalar eta12 = _q1.y() / _q.y();
        const Scalar eta21 = _q2.x() / _q.y();
        const Scalar eta22 = _q2.y() / _q.y();
        G_T << 0, 0, eta / ln, eta12 / 2, -eta11 / 2, 0, 0, 0, -eta / ln, eta22 / 2, -eta21 / 2, 0,
            0, 0, 1 / ln, 0, 0, 0, 0, 0, -1 / ln, 0, 0, 0,
            0, -1 / ln, 0, 0, 0, 0, 0, 1 / ln, 0, 0, 0, 0;
    }

    using Mat6_12 = Eigen::Matrix<Scalar, 6, 12>;
    Mat6_12 P{Mat6_12::Zero()};
    P.block<3, 3>(0, 3) = Mat3::Identity();
    P.block<3, 3>(3, 9) = Mat3::Identity();
    DEBUG_ONLY(cout << "P\n"
                    << P << endl;)
    P.block<3, 12>(0, 0) -= G_T;
    P.block<3, 12>(3, 0) -= G_T;

    DEBUG_ONLY(cout << "P\n"
                    << P << endl;)

    using Mat7_12 = Eigen::Matrix<Scalar, 7, 12>;
    Mat7_12 Ba{Mat7_12::Zero()};
    Ba.block<1, 12>(0, 0) = r_T;

    DEBUG_ONLY(cout << "Ba\n"
                    << Ba << endl;)

    Ba.block<6, 12>(1, 0) = P * E.transpose();

    DEBUG_ONLY(cout << "Ba\n"
                    << Ba << endl;)

    // Check this!!
    const Vec7 fl_o = calc_element_forces_local(ri_e, ro_e, l0, youngs, G, ul, theta_l1.x(), theta_l1.z(), -theta_l1.y(), theta_l2.x(), theta_l2.z(), -theta_l2.y());
    Vec7 fl{fl_o[3], fl_o[0], fl_o[1], fl_o[2], fl_o[4], fl_o[5], fl_o[6]};

    Mat7_12 B = Bl * Ba;
    DEBUG_ONLY(
        cout << "B\n"
             << B << endl;);
    // Vec12 dp;
    // dp << 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    // Vec7 dpl = B * dp;
    // cout << "dpl \n"
    //      << dpl << endl;

    const Vec12 fg = B.transpose() * fl;

    DEBUG_ONLY(cerr << "fl:\n"
                    << fl << endl;
               cerr << "fg:\n"
                    << fg << endl;)
    // R_int_trans[ie] += fg.segment(0, 3);
    // R_int_rot[ie] += fg.segment(3, 3);
    // R_int_trans[ie + 1] += fg.segment(6, 3);
    // R_int_rot[ie + 1] += fg.segment(9, 3);
    int a = 1;
}

inline void calc_element_inner_forces_crisfield(Index ie, const Vec3 *__restrict__ X, const Vec3 *__restrict__ d_trans,
                                                const Quaternion *__restrict__ d_rot, Vec3 *__restrict__ R_int_trans,
                                                Vec3 *__restrict__ R_int_rot, Scalar ri_e, Scalar ro_e, Scalar E, Scalar G)
{
    DEBUG_ONLY(cout << "!!!!!!!!!!!CRISFIELD!!!!!!!\n");
    DEBUG_ONLY(if (ie == 0) cout << "\n!!!!!!!!!!!! n = " << n_glob << " !!!!!!!!!!\n";);
    // if (n_glob == 2)
    // {
    //     int a = 1;
    // }

    // assert(i&e < X.size() - 1);
    const Vec3 &X1 = X[ie];
    const Vec3 &X2 = X[ie + 1];
    const Vec3 &d1 = d_trans[ie];
    const Vec3 &d2 = d_trans[ie + 1];

    DEBUG_ONLY(
        cout << "X1 " << X1.transpose() << endl;
        cout << "X2 " << X2.transpose() << endl;);

    DEBUG_ONLY(
        cout << "d1 " << d1.transpose() << endl;
        cout << "d2 " << d2.transpose() << endl;);
    const Mat3 T = d_rot[ie].to_matrix();
    const Vec3 &t1 = T.col(0);
    const Vec3 &t2 = T.col(1);
    const Vec3 &t3 = T.col(2);

    const Mat3 U = d_rot[ie + 1].to_matrix();
    const Vec3 &u1 = U.col(0);
    const Vec3 &u2 = U.col(1);
    const Vec3 &u3 = U.col(2);

    DEBUG_ONLY(
        cout << "U from quat\n"
             << U << endl;

        cout << "T from quat\n"
             << T << endl;);

    // calculate the first unit vector
    const Scalar l0 = (X2 - X1).norm();
    const Scalar ln = (X2 + d2 - (X1 + d1)).norm();

    const Vec3 e1 = (X2 + d2 - (X1 + d1)) / ln;
    assert(is_close(e1.norm(), 1.0));
    DEBUG_ONLY(
        cout << "e1 = " << e1.transpose() << endl;);
    // Calculate intermediate rotation between triads
    const Mat3 DeltaR = U * T.transpose();
    Quaternion qDeltaR;
    qDeltaR.from_matrix(DeltaR);
    assert(!is_close(qDeltaR.q0, 0.0, 0.1));        // if q0 is 0, gamma_half becomes singular
    const Vec3 gamma_half = qDeltaR.q / qDeltaR.q0; // eq 16.34 divided by 2

    const Mat3 S = skew_symmetric(gamma_half);
    const Mat3 DeltaR_m = Mat3::Identity() + 1 / (1 + 0.25 * gamma_half.dot(gamma_half)) * (S + 0.5 * S * S);
    assert(is_orthogonal(DeltaR_m));
    const Mat3 R_ = DeltaR_m * T;
    const Vec3 &r1 = R_.col(0);
    const Vec3 &r2 = R_.col(1);
    const Vec3 &r3 = R_.col(2);

    // computing element unit vectors e2 and e3 by rotating the unit vector
    // r1 on to e1.
    // const Vec3 e2 = r2 - r2.dot(e1) / (1 + r1.dot(e1)) * (e1 + r1);
    // const Vec3 e3 = r3 - r3.dot(e1) / (1 + r1.dot(e1)) * (e1 + r1);

    //    cout << "r1.dot.e1 " << e1.dot(r1) << endl;
    const Vec3 e2 = r2 - 0.5 * r2.dot(e1) * (e1 + r1);
    const Vec3 e3 = r3 - 0.5 * r3.dot(e1) * (e1 + r1);

    // const Vec3 e3 = e1.cross(e2);

    Scalar th = 2 * acos(qDeltaR.q0);
    Vec3 Om = 2 * qDeltaR.q / qDeltaR.q0;
    Vec3 Om_half = Om / 2;
    Scalar om_half = Om_half.norm();
    Scalar th_half = 2 * atan(om_half / 2);
    DEBUG_ONLY(
        cout << "th half " << th_half * 180 / M_PI << endl;
        cout << "q_DeltaR " << qDeltaR << endl;
        cout << "DeltaR_m \n"
             << DeltaR_m << endl;
        cout << "th " << th * 180 / M_PI << endl;
        cout << "R_\n"
             << R_ << endl;

        cout << "e2 = " << e2.transpose() << endl;
        cout << "e3 = " << e3.transpose() << endl;);
    // compute local displacements from global displacements.
    // Taking u_l = ln-l0 is not recommended since subtracting two large
    // numbers may be inaccurate with limited precision. Better to adopt
    // ul = ln - l0 = (ln - l0)*(ln + l0)/(ln + l0) = (ln^2 - l0^2)/(ln + l0)

    const Scalar ul = (ln * ln - l0 * l0) / (ln + l0);
    assert(is_close(ul, ln - l0));

    const Scalar theta_l1 = asin(0.5 * (-t3.dot(e2) + t2.dot(e3)));
    const Scalar theta_l2 = asin(0.5 * (-t2.dot(e1) + t1.dot(e2)));
    const Scalar theta_l3 = asin(0.5 * (-t3.dot(e1) + t1.dot(e3)));
    const Scalar theta_l4 = asin(0.5 * (-u3.dot(e2) + u2.dot(e3)));
    const Scalar theta_l5 = asin(0.5 * (-u2.dot(e1) + u1.dot(e2)));
    const Scalar theta_l6 = asin(0.5 * (-u3.dot(e1) + u1.dot(e3)));

    DEBUG_ONLY(
        cout << "Delta torsion [deg] " << 180 / M_PI * (theta_l4 - theta_l1) << endl;
        cout << "ul " << ul << endl;
        cout << "theta_l1 " << theta_l1 << endl;
        cout << "theta_l2 " << theta_l2 << endl;
        cout << "theta_l3 " << theta_l3 << endl;
        cout << "theta_l4 " << theta_l4 << endl;
        cout << "theta_l5 " << theta_l5 << endl;
        cout << "theta_l6 " << theta_l6 << endl;);
    // cout << "n_glob " << n_glob << endl;
    // assert(theta_l1 == 0);
    // assert(theta_l4 == 0);

#define MAX_ANGLE 85 * M_PI / 180
    assert(abs(theta_l1) < MAX_ANGLE);
    assert(abs(theta_l2) < MAX_ANGLE);
    assert(abs(theta_l3) < MAX_ANGLE);
    assert(abs(theta_l4) < MAX_ANGLE);
    assert(abs(theta_l5) < MAX_ANGLE);
    assert(abs(theta_l6) < MAX_ANGLE);
#undef MAX_ANGLE

    using Mat7_12 = Eigen::Matrix<Scalar, 12, 7>;
    Mat7_12 F_transpose{Mat7_12::Zero()}; // Transpose of F where zero rows in F are excluded
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

    Eigen::Matrix<Scalar, 12, 3> Lr2{}, Lr3{};
    Lr2 << L1r2, L2r2, -L1r2, L2r2; // Ikke transponer?
    Lr3 << L1r3, L2r3, -L1r3, L2r3;

    Vec12 h1{Vec12::Zero()}, h2{Vec12::Zero()}, h3{Vec12::Zero()}, h4{Vec12::Zero()}, h5{Vec12::Zero()}, h6{Vec12::Zero()};
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
    DEBUG_ONLY(
        cout << "Crisfield check\n";);
    Vec12 dp;
    dp << 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    Vec7 dpl = F_transpose.transpose() * dp;
    DEBUG_ONLY(
        cout << "dpl cf\n"
             << dpl << endl;);
    // cout << "F^T\n"
    //      << F.transpose() << endl;

    /*Calculate local internal forces based on linear 3D beam theory*/
    const Vec7 R_int_e_l = calc_element_forces_local(ri_e, ro_e, l0, E, G, ul,
                                                     theta_l1, theta_l2, theta_l3, theta_l4, theta_l5, theta_l6);

    DEBUG_ONLY(
        cout << "R_int_e_l:\n"
             << R_int_e_l << endl;);

    const Vec12 R_int_e = F_transpose * R_int_e_l;
    assert(R_int_e.allFinite());
    DEBUG_ONLY(
        cout << "R_int_e:\n"
             << R_int_e << endl;);

    // R_int_trans[ie] += R_int_e.segment(0, 3);
    // R_int_rot[ie] += R_int_e.segment(3, 3);
    // R_int_trans[ie + 1] += R_int_e.segment(6, 3);
    // R_int_rot[ie + 1] += R_int_e.segment(9, 3);
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
    Scalar M1 = G * K * (theta_1l - theta_4l) / l0;
    Scalar M4 = -M1;

    // M1*=-1;
    // M4*=-1;

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
    DEBUG_ONLY(
        print_std_vector(beam_system.R_int_trans, "R int trans");
        print_std_vector(beam_system.R_int_rot, "R int rot"););

#pragma omp parallel for
    /*even elements*/
    for (Index ie = 0; ie < Ne; ie += 2)
    {
        // if (ie % 2 == 0)
        // {
        calc_element_inner_forces(ie, X.data(), beam_system.d_trans.data(), beam_system.d_rot.data(),
                                  beam_system.R_int_trans.data(), beam_system.R_int_rot.data(),
                                  geometry.ri_e(ie), geometry.ro_e(ie), E, G);
        calc_element_inner_forces_battini(ie, X.data(), beam_system.d_trans.data(), beam_system.d_rot.data(),
                                          beam_system.R_int_trans.data(), beam_system.R_int_rot.data(),
                                          geometry.ri_e(ie), geometry.ro_e(ie), E, G);
        // }
    }
#pragma omp parallel for
    /*Odd elements*/
    for (Index ie = 1; ie < Ne; ie += 2)
    {

        // if (ie % 2 == 1)
        // {
        calc_element_inner_forces(ie, X.data(), beam_system.d_trans.data(), beam_system.d_rot.data(),
                                  beam_system.R_int_trans.data(), beam_system.R_int_rot.data(),
                                  geometry.ri_e(ie), geometry.ro_e(ie), E, G);

        calc_element_inner_forces_battini(ie, X.data(), beam_system.d_trans.data(), beam_system.d_rot.data(),
                                          beam_system.R_int_trans.data(), beam_system.R_int_rot.data(),
                                          geometry.ri_e(ie), geometry.ro_e(ie), E, G);
        // }
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

// inline void velocity_update_partial_OLD(Scalar dt, Index N, const Scalar *__restrict__ M, const Vec3 *__restrict__ J_u,
//                                         const Vec3 *__restrict__ R_int_trans, const Vec3 *__restrict__ R_int_rot,
//                                         const Vec3 *__restrict__ R_ext_trans, const Vec3 *__restrict__ R_ext_rot,
//                                         Vec3 *__restrict__ v_trans, Vec3 *__restrict__ v_rot)
// {
// #pragma omp parallel
//     {
// #pragma omp for nowait
//         for (Index i = 0; i < N; i++)
//         {
//             v_trans[i] += 0.5 * dt * (R_ext_trans[i] - R_int_trans[i]) / M[i];

//             DEBUG_ONLY(
//                 if (i == 1) cout << "v_trans:\n"
//                                  << v_trans[i] << endl;);
//         }
// #pragma omp for
//         for (Index i = 0; i < N; i++)
//         {
//             /*omega_u_dot = J_u^-1 * (R_rot_u - S(omega_u)*J_u*omega_u)*/
//             const Vec3 R_rot_u = R_ext_rot[i] - R_int_rot[i];

//             DEBUG_ONLY(if (i == 1) cout
//                        << "R_ext_rot\n"
//                        << R_ext_rot[i] << endl
//                        << "R_int_rot\n"
//                        << R_int_rot[i] << endl
//                        << "R_rot\n"
//                        << R_rot_u << endl);
//             Vec3 &omega_u = v_rot[i];

//             omega_u = {1.424, -2.13, 3.11331};

//             Mat3 JJu = J_u[i].asDiagonal();
//             Vec3 Ju = J_u[i];
//             Vec3 rot_term_orig = omega_u.cross(Vec3{J_u[i].array() * omega_u.array()});

//             Vec3 rot_term_new = skew_symmetric(omega_u) * JJu * omega_u;

//             cout << "rot_term_orig " << rot_term_orig << endl;
//             cout << "rot_term_new " << rot_term_new << endl;

//             Vec3 omega_u_dot_new;
//             omega_u_dot_new.x() = (R_rot_u.x() - (Ju.z() - Ju.y()) * omega_u.y() * omega_u.z()) / Ju.x();
//             omega_u_dot_new.y() = (R_rot_u.y() - (Ju.x() - Ju.z()) * omega_u.x() * omega_u.z()) / Ju.y();
//             omega_u_dot_new.z() = (R_rot_u.z() - (Ju.y() - Ju.x()) * omega_u.x() * omega_u.y()) / Ju.z();

//             const Vec3 omega_u_dot = (R_rot_u - rot_term_orig).array() / J_u[i].array();

//             cout << "omega_u_dot\n"
//                  << omega_u_dot << endl;
//             cout << "omega_u_dot_new\n"
//                  << omega_u_dot_new << endl;
//             // cout << "diff: \n"
//             //      << omega_u_dot - omega_u_dot_new << endl;

//             // const Vec3 omega_u_dot = (R_rot_u - omega_u.cross(Vec3{J_u[i].array() * omega_u.array()})).array() / J_u[i].array();
//             //  Scalar tol = 0.00001;
//             //  if (i != 0 && abs(omega_u[0]) > tol)
//             //  {
//             //      cout << "omega_u \n"
//             //           << omega_u << endl;
//             //      assert(false);
//             //  }
//             omega_u += 0.5 * dt * omega_u_dot;
//             DEBUG_ONLY(
//                 if (i == 1) cout << "omega_u:\n"
//                                  << v_rot[i] << endl;);
//             // if (n_glob == 10000 && i > 10)
//             //     omega_u[0] = 0.1;
//         }
//     }
// }

inline void displacement_update(Scalar dt, Index N, Vec3 *__restrict__ v_trans, Vec3 *__restrict__ v_rot, Vec3 *__restrict__ d_trans,
                                Quaternion *__restrict__ d_rot)
{
#pragma omp parallel
    {
#pragma omp for nowait
        for (Index i = 0; i < N; i++)
        {
            d_trans[i] += dt * v_trans[i];
            DEBUG_ONLY(
                if (i == 1) cout << "d_trans:\n"
                                 << d_trans[i] << endl;);
        }
#pragma omp for
        for (Index i = 0; i < N; i++)
        {
            // Mat3 U = d_rot[i].to_matrix();
            // U = U * (Mat3::Identity() - 0.5 * dt * skew_symmetric(v_rot[i])).inverse() *
            //     (Mat3::Identity() + 0.5 * dt * skew_symmetric(v_rot[i]));
            // assert(U.allFinite());
            // assert(is_orthogonal(U));
            // d_rot[i].from_matrix(U);
            // assert(is_close(d_rot[i].norm(), 1.0));

            Quaternion &q = d_rot[i];
            const Vec3 &omega_u = v_rot[i];

            const Vec3 omega = q.rotate_vector(omega_u); // perhaps possible to simplify this and not having to first convert omega to the global frame

            // const Vec3 omega = q.rotate_vector_reversed(omega_u); // perhaps possible to simplify this and not having to first convert omega to the global frame
            if (i == 1)
            {
                cout << "i " << i << endl;
                cout << "omega " << omega << endl;
                cout << "omega mag " << omega.norm() << endl;
            }
            // assert(is_close(omega_u.y(), 0, 1e-4) && is_close(omega_u.z(), 0, 1e-4));
            q.compound_rotate(dt * omega);
            Scalar norm = q.norm();
            assert(is_close(norm, 1.0));
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
        if (i == 1)
        {
            cout << "R_ext_rot before\n"
                 << R_ext_rot[i] << endl;
        }
        //     Mat3 U = d_rot[i].to_matrix();
        //     R_int[i].rot = U.transpose() * R_int[i].rot;
        //     R_ext[i].rot = U.transpose() * R_ext[i].rot;
        const Quaternion &q = d_rot[i];

        R_int_rot[i] = q.rotate_vector_reversed(R_int_rot[i]);
        R_ext_rot[i] = q.rotate_vector_reversed(R_ext_rot[i]);

        if (i == 1)
        {
            cout << "R_ext_rot after\n"
                 << R_ext_rot[i] << endl;
        }
    }
}

inline void add_mass_proportional_rayleigh_damping(Index N, Scalar alpha, const Scalar *__restrict__ M,
                                                   const Vec3 *__restrict__ v_trans, Vec3 *__restrict__ R_int_trans,
                                                   const Vec3 *__restrict__ J_u, const Quaternion *__restrict__ d_rot,
                                                   const Vec3 *__restrict__ v_rot, Vec3 *__restrict__ R_int_rot)
{
#pragma omp parallel for
    for (Index i = 0; i < N; i++)
    {
        R_int_trans[i] += alpha * M[i] * v_trans[i];
    }
// Test with rotational dofs also?
#pragma omp parallel for
    for (Index i = 0; i < N; i++)
    {
        Scalar alpha_rot = 1 * alpha;
        Vec3 R_damp_rot = alpha_rot * J_u[i].array() * v_rot[i].array();
        R_damp_rot = d_rot->rotate_vector(R_damp_rot);
        // if (i == 10 && n_glob % 100 == 0)
        // {
        //     cout << "M_w_1 prior: " << R_int_rot[i].x() << endl;
        //     cout << "M damp w_1: " << R_damp_rot.x() << endl;
        // }
        R_int_rot[i] += R_damp_rot;
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

inline int is_node_within_hole_segment(Index i, const Vec3 &x_hole_A, const Vec3 &x_hole_B,
                                       const Vec3 &X, const Vec3 &d_trans)
{
    const Vec3 t = (x_hole_B - x_hole_A).normalized();
    const Vec3 x = X + d_trans;
    const Scalar l_hole_seg = (x_hole_B - x_hole_A).norm();
    const Scalar dist = (x - x_hole_A).dot(t);
    const bool is_between = dist >= -SMALL_SCALAR && dist <= l_hole_seg + SMALL_SCALAR;
    if (is_between)
    {
        return 0;
    }
    else if (dist < 0)
    {
        return -1;
    }
    else
    {
        assert(dist > 0);
        return 1;
    }
}

inline void update_hole_contact_indices(const Index N, const Vec3 *__restrict__ x_hole,
                                        Index *__restrict__ hole_index, const Vec3 *__restrict__ X,
                                        const Vec3 *__restrict__ d_trans)
{
    const Index Ne = N - 1;
    /*Update hole indices*/
    for (Index i = 0; i < N; i++)
    {
        int hi = hole_index[i];
        int is_between = is_node_within_hole_segment(i, x_hole[hi], x_hole[hi + 1], X[i], d_trans[i]);
        if (is_between == 0)
        {
            continue;
        }
        else if (is_between == -1)
        {
            /*Search backwards*/
            hi--;
            assert(hi >= 0);
            while (is_node_within_hole_segment(i, x_hole[hi], x_hole[hi + 1], X[i], d_trans[i]) != 0)
            {
                assert(is_node_within_hole_segment(i, x_hole[hi], x_hole[hi + 1], X[i], d_trans[i]) == -1);
                hi--;
                assert(hi >= 0);
            }
            hole_index[i] = hi;
        }
        else
        {
            assert(is_between == 1);
            /*Search forwards*/
            hi++;
            assert(hi < Ne);
            while (is_node_within_hole_segment(i, x_hole[hi], x_hole[hi + 1], X[i], d_trans[i]) != 0)
            {
                assert(is_node_within_hole_segment(i, x_hole[hi], x_hole[hi + 1], X[i], d_trans[i]) == 1);
                hi++;
                assert(hi < Ne);
            }
            hole_index[i] = hi;
        }
    }
}

inline void calc_contact_forces(const Config &config, const Index N, const Vec3 *__restrict__ x_hole,
                                Index *__restrict__ hole_index, const Scalar *__restrict__ r_hole,
                                const Scalar *__restrict__ r_outer_string, const Vec3 *__restrict__ X,
                                const Vec3 *__restrict__ d_trans, const Quaternion *__restrict__ d_rot,
                                Vec3 *__restrict__ R_ext)
{
    update_hole_contact_indices(N, x_hole, hole_index, X, d_trans);

    const Index Ne = N - 1;
    for (Index i = 0; i < N; i++)
    {
        const Index hi = hole_index[i];
        assert(is_node_within_hole_segment(i, x_hole[hi], x_hole[hi + 1], X[i], d_trans[i]) == 0);
        assert(hi < Ne);
        const Vec3 &x = X[i] + d_trans[i];
        const Vec3 &x_hole_A = x_hole[hi];
        const Vec3 &x_hole_B = x_hole[hi + 1];
        const Vec3 t = (x_hole_B - x_hole_A).normalized();
        const Vec3 x_center = x_hole_A + (x - x_hole_A).dot(t) * t;
        const Scalar d_center = (x - x_center).norm();
        const Scalar delta = d_center + r_outer_string[i] - r_hole[i];
        if (delta <= 0.0)
        {
            continue;
        }
        const Vec3 n = (x - x_center).normalized();
    }
}

inline void set_simple_bc(const Config &config, const Geometry &geometry, BeamSystem &beam_system)
{
    Index N = geometry.get_N();

    vector<Vec3> d_trans = beam_system.d_trans;
    vector<Quaternion> &d_rot = beam_system.d_rot;
    vector<Vec3> &v_trans = beam_system.v_trans;
    vector<Vec3> &v_rot = beam_system.v_rot;

    assert(N == d_trans.size() && N == d_rot.size() && N == v_trans.size() && N == v_rot.size());

    switch (config.bc_case)
    {
    case BC_Case::NONE:
        return;
    case BC_Case::CANTILEVER:
        d_trans[0] = {0, 0, 0};
        v_trans[0] = {0, 0, 0};
        // d_rot[0].from_matrix();
        d_rot[0] = config.bc_orientation_base;
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