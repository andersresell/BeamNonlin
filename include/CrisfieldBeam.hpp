#pragma once

#include "Config.hpp"
#include "SolverUtils.hpp"

namespace CrisfieldBeam {

Vec12 calc_element_inner_forces(const Index ie, const vector<Vec3> &X, const vector<Vec3> &d_trans,
                                const vector<Quaternion> &d_rot, const Scalar youngs, const Scalar G, const Scalar I_2,
                                const Scalar I_3, const Scalar A, const Scalar J, const Scalar beta_rayleigh,
                                const vector<Vec3> &v_trans, const vector<Vec3> &v_rot);

Vec7 calc_element_forces_local(const Scalar l0, const Scalar youngs, const Scalar G, const Scalar I_2, const Scalar I_3,
                               const Scalar A, const Scalar J, const Vec7 &d_l);

inline void calc_element_kinematics(const Index ie, const vector<Vec3> &X, const vector<Vec3> &d_trans,
                                    const vector<Quaternion> &d_rot, Scalar &l0, Scalar &ln, Mat3 &T, Mat3 &U, Mat3 &E,
                                    Mat3 &R_) {
    const Vec3 &X1 = X[ie];
    const Vec3 &X2 = X[ie + 1];
    const Vec3 &d1 = d_trans[ie];
    const Vec3 &d2 = d_trans[ie + 1];

    T = d_rot[ie].to_matrix();
    U = d_rot[ie + 1].to_matrix();

    l0 = (X2 - X1).norm();
    ln = (X2 + d2 - (X1 + d1)).norm();

    // Calculate intermediate rotation between triads
    const Mat3 DeltaR = U * T.transpose();
    Quaternion qDeltaR;
    qDeltaR.from_matrix(DeltaR);
    assert(!is_close(qDeltaR.q0, 0.0, 0.1));        // if q0 is 0, gamma_half becomes singular
    const Vec3 gamma_half = qDeltaR.q / qDeltaR.q0; // eq 16.34 divided by 2

    const Mat3 S = skew(gamma_half);
    const Mat3 DeltaR_m = Mat3::Identity() + 1.0 / (1.0 + 0.25 * gamma_half.dot(gamma_half)) * (S + 0.5 * S * S);
    assert(is_orthogonal(DeltaR_m));
    R_ = DeltaR_m * T;
    const Vec3 &r1 = R_.col(0);
    const Vec3 &r2 = R_.col(1);
    const Vec3 &r3 = R_.col(2);

    E.col(0) = (X2 + d2 - (X1 + d1)) / ln;
    const Vec3 &e1 = E.col(0);
    assert(is_close(e1.norm(), 1.0));
    // computing element unit vectors e2 and e3 by rotating the unit vector
    // r1 on to e1.
    E.col(1) = r2 - r2.dot(e1) / (1 + r1.dot(e1)) * (e1 + r1);
    E.col(2) = r3 - r3.dot(e1) / (1 + r1.dot(e1)) * (e1 + r1);
    assert(is_orthogonal(E));
}

inline Vec7 calc_disp_local(const Scalar l0, const Scalar ln, const Mat3 &T, const Mat3 &U, const Mat3 &E) {
    const Vec3 &t1 = T.col(0);
    const Vec3 &t2 = T.col(1);
    const Vec3 &t3 = T.col(2);

    const Vec3 &u1 = U.col(0);
    const Vec3 &u2 = U.col(1);
    const Vec3 &u3 = U.col(2);

    const Vec3 &e1 = E.col(0);
    const Vec3 &e2 = E.col(1);
    const Vec3 &e3 = E.col(2);

    // compute local displacements from global displacements.
    // Taking u_l = ln-l0 is not recommended since subtracting two large
    // numbers may be inaccurate with limited precision. Better to adopt
    // ul = ln - l0 = (ln - l0)*(ln + l0)/(ln + l0) = (ln^2 - l0^2)/(ln + l0)
    const Vec7 d_l = {
        (ln * ln - l0 * l0) / (ln + l0),        asin(0.5 * (-t3.dot(e2) + t2.dot(e3))),
        asin(0.5 * (-t2.dot(e1) + t1.dot(e2))), asin(0.5 * (-t3.dot(e1) + t1.dot(e3))),
        asin(0.5 * (-u3.dot(e2) + u2.dot(e3))), asin(0.5 * (-u2.dot(e1) + u1.dot(e2))),
        asin(0.5 * (-u3.dot(e1) + u1.dot(e3))),
    };
    return d_l;
}

using Mat7_12 = Eigen::Matrix<Scalar, 12, 7>;
inline Mat7_12 calc_F_transpose(const Mat3 &T, const Mat3 &U, const Mat3 &E, const Mat3 &R_, const Vec7 &d_l,
                                const Scalar ln) {
    constexpr Index f4 = 0;
    constexpr Index f5 = 1;
    constexpr Index f6 = 2;
    constexpr Index f7 = 3;
    constexpr Index f10 = 4;
    constexpr Index f11 = 5;
    constexpr Index f12 = 6;

    const Vec3 &t1 = T.col(0);
    const Vec3 &t2 = T.col(1);
    const Vec3 &t3 = T.col(2);

    const Vec3 &u1 = U.col(0);
    const Vec3 &u2 = U.col(1);
    const Vec3 &u3 = U.col(2);

    const Vec3 &e1 = E.col(0);
    const Vec3 &e2 = E.col(1);
    const Vec3 &e3 = E.col(2);

    const Vec3 &r1 = R_.col(0);
    const Vec3 &r2 = R_.col(1);
    const Vec3 &r3 = R_.col(2);

    Mat7_12 F_transpose; // Transpose of F where zero rows in F are excluded

    // calculate F connecting infinitesimal global and local variable
    // optimize zero rows later

    F_transpose.col(f7) << -e1, Vec3::Zero(), e1, Vec3::Zero(); //(17.19)

    const Mat3 A = 1.0 / ln * (Mat3::Identity() - e1 * e1.transpose());
    assert(is_close((A - A.transpose()).norm(), 0.0));

    const Mat3 L1r2 = 0.5 * r2.dot(e1) * A + 0.5 * A * r2 * (e1 + r1).transpose();
    const Mat3 L1r3 = 0.5 * r3.dot(e1) * A + 0.5 * A * r3 * (e1 + r1).transpose();
    const Mat3 L2r2 =
        0.5 * skew(r2) - 0.25 * r2.transpose() * e1 * skew(r1) - 0.25 * skew(r2) * e1 * (e1 + r1).transpose();
    const Mat3 L2r3 =
        0.5 * skew(r3) - 0.25 * r3.transpose() * e1 * skew(r1) - 0.25 * skew(r3) * e1 * (e1 + r1).transpose();

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

    F_transpose.col(f4) = 1 / (2 * cos(d_l[1])) * (Lr3 * t2 - Lr2 * t3 + h1);
    F_transpose.col(f5) = 1 / (2 * cos(d_l[2])) * (Lr2 * t1 + h2);
    F_transpose.col(f6) = 1 / (2 * cos(d_l[3])) * (Lr3 * t1 + h3);
    F_transpose.col(f10) = 1 / (2 * cos(d_l[4])) * (Lr3 * u2 - Lr2 * u3 + h4);
    F_transpose.col(f11) =
        1 / (2 * cos(d_l[5])) * (Lr2 * u1 + h5); // seems to be an error in the book for f11 and f12. it should be +
                                                 // in front of h6 and h6 (it's + in the paper)
    F_transpose.col(f12) = 1 / (2 * cos(d_l[6])) * (Lr3 * u1 + h6);
    return F_transpose;
}

} // namespace CrisfieldBeam
