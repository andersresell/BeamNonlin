#pragma once
#include "Utils.hpp"
#include "Containers.hpp"

inline void calc_element_contribution(Index ie, const vector<Vec3> &nodes, vector<GlobalDisp> &disp_g,
                                      vector<Scalar> &R_int)
{
    assert(ie < nodes.size() - 1);
    const Vec3 &x1 = nodes[ie];
    const Vec3 &x2 = nodes[ie + 1];
    const Vec3 &u1 = disp_g[ie].u_trans;
    const Vec3 &u2 = disp_g[ie + 1].u_trans;

    const Mat3 T = disp_g[ie].u_quat.to_triad();
    const Vec3 t1 = T.col(0); // consider optimizing by using T directly
    const Vec3 t2 = T.col(1);
    const Vec3 t3 = T.col(2);

    const Mat3 U = disp_g[ie + 1].u_quat.to_triad();
    const Vec3 u1 = U.col(0); // consider optimizing by using U directly
    const Vec3 u2 = U.col(1);
    const Vec3 u3 = U.col(2);

    // calculate the first unit vector
    const Scalar l0 = (x2 - x1).norm();
    const Scalar ln = (x2 + u2 - (x1 + u1)).norm();
    const Vec3 e1 = (x2 + u2 - (x1 + u1)) / ln;
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

    // calculate F connecting infinitesimal global and local variables
}

inline void assemble(const Config &config, const Geometry &geometry)
{

    Index Ne = geometry.get_Ne();

    for (Index ie = 0; ie < Ne; ie++)
    {
    }
}