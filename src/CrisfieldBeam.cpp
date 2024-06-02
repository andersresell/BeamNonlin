#include "../include/CrisfieldBeam.hpp"

namespace CrisfieldBeam {
static Vec12 calc_approx_rayleigh_beta_damping(const Scalar beta, const Scalar A, const Scalar I_2, const Scalar I_3,
                                               const Scalar J, const Scalar l0, const Scalar youngs, const Scalar G,
                                               const Mat3 &E, const Mat3 &U1, const Mat3 &U2, const Vec3 &v1,
                                               const Vec3 &v2, const Vec3 &omega1_u, const Vec3 &omega2_u);

Vec12 calc_element_inner_forces(const Index ie, const vector<Vec3> &X, const vector<Vec3> &d_trans,
                                const vector<Quaternion> &d_rot, const Scalar youngs, const Scalar G, const Scalar I_2,
                                const Scalar I_3, const Scalar A, const Scalar J, const Scalar beta_rayleigh,
                                const vector<Vec3> &v_trans, const vector<Vec3> &v_rot) {

    Scalar l0, ln;
    Mat3 T, U, E, R_;
    calc_element_kinematics(ie, X, d_trans, d_rot, l0, ln, T, U, E, R_);

    const Vec7 d_l = calc_disp_local(l0, ln, T, U, E);

    const Vec7 R_int_e_l = calc_element_forces_local(l0, youngs, G, I_2, I_3, A, J, d_l);

    const Mat7_12 F_transpose = calc_F_transpose(T, U, E, R_, d_l, ln);

    Vec12 R_int_e = F_transpose * R_int_e_l;
    assert(R_int_e.allFinite());

    if (beta_rayleigh > 0) {
        const Vec12 R_damp = calc_approx_rayleigh_beta_damping(beta_rayleigh, A, I_2, I_3, J, l0, youngs, G, E, U, T,
                                                               v_trans[ie], v_trans[ie + 1], v_rot[ie], v_rot[ie + 1]);
        DEBUG_ONLY(cout << "R_damp\n" << R_damp << endl;);
        R_int_e += R_damp;
    } else {
        assert(beta_rayleigh == 0.0);
    }

    return R_int_e;
}

Vec7 calc_element_forces_local(const Scalar l0, const Scalar youngs, const Scalar G, const Scalar I_2, const Scalar I_3,
                               const Scalar A, const Scalar J, const Vec7 &d_l) {

    /*Normal force (F1)
     [[F1], = A*E/l0[[ 1 -1],*[[0],
      [F4]]          [-1  1]]  [ul]]
    */
    const Scalar ul = d_l[0];
    const Scalar theta_1l = d_l[1];
    const Scalar theta_2l = d_l[2];
    const Scalar theta_3l = d_l[3];
    const Scalar theta_4l = d_l[4];
    const Scalar theta_5l = d_l[5];
    const Scalar theta_6l = d_l[6];

    const Scalar F1 = A * youngs * (-ul) / l0;
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
    const Scalar M1 = G * J * (theta_1l - theta_4l) / l0;
    const Scalar M4 = -M1;

    /*Bending: Euler bernoulli with only angle dofs:
    k = EI/L [[4 2],
              [2 4]]
    */
    const Scalar M2 = youngs * I_2 / l0 * (4 * theta_2l + 2 * theta_5l);
    const Scalar M5 = youngs * I_2 / l0 * (2 * theta_2l + 4 * theta_5l);

    const Scalar M3 = youngs * I_3 / l0 * (4 * theta_3l + 2 * theta_6l);
    const Scalar M6 = youngs * I_3 / l0 * (2 * theta_3l + 4 * theta_6l);

    // Vec12 R_int_l = {F1, 0, 0, M1, M2, M3, F4, 0, 0, M4, M5, M6};
    const Vec7 R_int_e_l = {M1, M2, M3, F4, M4, M5, M6};
    // cout << "fix\n";
    // //Vec12 R_int_l = {0, 0, 0, 0, M2, M3, 0, 0, 0, 0, M5, M6};
    // Vec12 R_int_l = {0, 0, 0, 0, M2, M3, 0, 0, 0, 0, M5, M6};
    assert(R_int_e_l.allFinite());
    return R_int_e_l;
}

static Vec12 calc_approx_rayleigh_beta_damping(
    const Scalar beta, const Scalar A, const Scalar I_2, const Scalar I_3, const Scalar J, const Scalar l0,
    const Scalar youngs, const Scalar G, const Mat3 &E, const Mat3 &U1, const Mat3 &U2, const Vec3 &v1, const Vec3 &v2,
    const Vec3 &omega1_u, const Vec3 &omega2_u) { /*Basing this on the formula R_damp_l = beta * k_l * v_l,
                                                 First the local velocity v_l is computed by rotating the global
                                                 components, then the local damping is computed, and finally the
                                                 resulting force is rotated to global components*/

    const Vec3 v1_l = E.transpose() * v1;
    const Vec3 v2_l = E.transpose() * v2;

    const Vec3 omega1_l = E.transpose() * U1 * omega1_u;
    const Vec3 omega2_l = E.transpose() * U2 * omega2_u;

    Mat4 k_EB = Mat4{{12, 6 * l0, -12, 6 * l0},
                     {6 * l0, 4 * l0 * l0, -6 * l0, 2 * l0 * l0},
                     {-12, -6 * l0, 12, -6 * l0},
                     {6 * l0, 2 * l0 * l0, -6 * l0, 4 * l0 * l0}} *
                youngs * I_2 / (l0 * l0 * l0);

    const Scalar K = J; // cicular cross-section
    const Mat2 k_tor = Mat2{{1, -1}, {-1, 1}} * G * K / l0;
    const Mat2 k_x = Mat2{{1, -1}, {-1, 1}} * youngs * A / l0;

    const Vec2 v_x = {v1_l.x(), v2_l.x()};
    const Vec2 v_tor = {omega1_l.x(), omega2_l.x()};
    const Vec4 v_y = {v1_l.y(), omega1_l.z(), v2_l.y(), omega2_l.z()};
    const Vec4 v_z = {v1_l.z(), -omega1_l.y(), v2_l.z(), -omega2_l.y()};

    const Vec2 f_x = beta * k_x * v_x;
    const Vec2 f_tor = beta * k_tor * v_tor;
    const Vec4 f_y = beta * k_EB * v_y;
    k_EB *= (I_3 / I_2);
    const Vec4 f_z = beta * k_EB * v_z;

    const Vec3 f1_l = {f_x[0], f_y[0], f_z[0]};
    const Vec3 m1_l = {f_tor[0], -f_z[1], f_y[1]};
    const Vec3 f2_l = {f_x[1], f_y[2], f_z[2]};
    const Vec3 m2_l = {f_tor[1], -f_z[3], f_y[3]};

    Vec12 R_damp;
    R_damp << E * f1_l, E * m1_l, E * f2_l, E * m2_l;
    return R_damp;
}
} // namespace CrisfieldBeam