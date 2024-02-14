#pragma once
#include "Utils.hpp"
#include "Containers.hpp"

inline void calc_element_contribution(Index ie, const vector<Vec3> &X, vector<GlobalDisp> &disp_g,
                                      vector<Vec6> &R_int)
{
    assert(ie < X.size() - 1);
    const Vec3 &X1 = X[ie];
    const Vec3 &X2 = X[ie + 1];
    const Vec3 &u1 = disp_g[ie].u_trans;
    const Vec3 &u2 = disp_g[ie + 1].u_trans;

    const Mat3 T = disp_g[ie].u_quat.to_triad();
    const Vec3 t1 = T.col(0);
    const Vec3 t2 = T.col(1);
    const Vec3 t3 = T.col(2);

    const Mat3 U = disp_g[ie + 1].u_quat.to_triad();
    const Vec3 u1 = U.col(0);
    const Vec3 u2 = U.col(1);
    const Vec3 u3 = U.col(2);

    // calculate the first unit vector
    const Scalar l0 = (X2 - X1).norm();
    const Scalar ln = (X2 + u2 - (X1 + u1)).norm();
    const Vec3 e1 = (X2 + u2 - (X1 + u1)) / ln;
    assert(is_close(e1.norm(), 1.0));

    // Calculate intermediate rotation between triads
    const Mat3 DeltaR = U * T.transpose();
    Quaternion qDeltaR;
    qDeltaR.from_triad(DeltaR);
    const Vec3 gamma_half = Vec3{qDeltaR.q1, qDeltaR.q2, qDeltaR.q3} / qDeltaR.q0; // eq 16.34 divided by 2
    const Mat3 S = skew_symmetric(gamma_half);
    const Mat3 DeltaR_m = Mat3::Identity() + 1 / (1 + 0.25 * gamma_half.dot(gamma_half)) * (S + 0.5 * S * S);
    const Mat3 R_ = DeltaR_m * T;
    const Vec3 r1 = R_.col(0);
    const Vec3 r2 = R_.col(1);
    const Vec3 r3 = R_.col(2);

    // computing element unit vectors e2 and e3 by rotating the unit vector
    // r1 on to e1.
    const Vec3 e2 = r2 - r2.dot(e1) / (1 + r1.dot(e1)) * (e1 + r1);
    const Vec3 e3 = r3 - r3.dot(e1) / (1 + r1.dot(e1)) * (e1 + r1);

    // compute local displacements from global displacements
    const Scalar ul = ln - l0;
    const Scalar theta_l1 = 0.5 * asin(-t3.dot(e2) + t2.dot(e3));
    const Scalar theta_l2 = 0.5 * asin(-t2.dot(e1) + t1.dot(e2));
    const Scalar theta_l3 = 0.5 * asin(-t3.dot(e1) + t1.dot(e3));
    const Scalar theta_l4 = 0.5 * asin(-u3.dot(e2) + u2.dot(e3));
    const Scalar theta_l5 = 0.5 * asin(-u2.dot(e1) + u1.dot(e2));
    const Scalar theta_l6 = 0.5 * asin(-u3.dot(e1) + u1.dot(e3));

    // calculate F connecting infinitesimal global and local variable
    // optimize zero rows later
    using Vec12 = Eigen::Vector<Scalar, 12>;
    Vec12 f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12;
    f1.setConstant(0);
    f2.setConstant(0);
    f3.setConstant(0);
    f8.setConstant(0);
    f9.setConstant(0);

    f7 << -e1, Vec3::Zero(), e1, Vec3::Zero(); //(17.19)

    const Vec3 A = 1 / ln * (Mat3::Identity() - e1 * e1.transpose());

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

    f4 = 1 / (2 * cos(theta_l1)) * (Lr3 * t2 - Lr2 * t3 + h1);
    f5 = 1 / (2 * cos(theta_l2)) * (Lr2 * t1 + h2);
    f6 = 1 / (2 * cos(theta_l3)) * (Lr3 * t1 + h3);
    f10 = 1 / (2 * cos(theta_l4)) * (Lr3 * u2 - Lr2 * u3 + h4);
    f11 = 1 / (2 * cos(theta_l5)) * (Lr2 * u1 - h5);
    f12 = 1 / (2 * cos(theta_l6)) * (Lr3 * u1 - h6);

    Eigen::Matrix<Scalar, 12, 12> F;
    F << f1.transpose(),
        f2.transpose(),
        f3.transpose(),
        f4.transpose(),
        f5.transpose(),
        f6.transpose(),
        f7.transpose(),
        f8.transpose(),
        f9.transpose(),
        f10.transpose(),
        f11.transpose(),
        f12.transpose();

    Vec12 R_int_l; // calc this using euler bernoulli e.l
    Vec12 R_int_e = F.transpose() * R_int_l;
    R_int[ie] += R_int_e.segment(0, 6);
    R_int[ie + 1] += R_int_e.segment(6, 6);
}

inline void assemble(const Config &config, const Geometry &geometry)
{

    Index Ne = geometry.get_Ne();

    for (Index ie = 0; ie < Ne; ie++)
    {
        calc_element_contribution(ie, )
    }
}