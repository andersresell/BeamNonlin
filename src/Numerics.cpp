

#include "../include/Numerics.hpp"

static void calc_forces(const Config &config, const Geometry &geometry, const Borehole &borehole, BeamSystem &beam);

template <Index i_first>
static void calc_inner_forces(const Config &config, const Geometry &geometry, BeamSystem &beam);

static void calc_element_inner_forces(const Index ie, const vector<Vec3> &X, const vector<Vec3> &d_trans,
                                      const vector<Quaternion> &d_rot, vector<Vec3> &R_int_trans,
                                      vector<Vec3> &R_int_rot, const Scalar youngs, const Scalar G, const Scalar I_2,
                                      const Scalar I_3, const Scalar A, const Scalar J, const Scalar beta_rayleigh,
                                      const vector<Vec3> &v_trans, const vector<Vec3> &v_rot);
static Vec7 calc_element_forces_local(const Scalar l0, const Scalar youngs, const Scalar G, const Scalar I_2,
                                      const Scalar I_3, const Scalar A, const Scalar J, const Scalar ul,
                                      const Scalar theta_1l, const Scalar theta_2l, const Scalar theta_3l,
                                      const Scalar theta_4l, const Scalar theta_5l, const Scalar theta_6l);
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
static Vec12 calc_approx_rayleigh_beta_damping(const Scalar beta, const Scalar A, const Scalar I_2, const Scalar I_3,
                                               const Scalar J, const Scalar l0, const Scalar youngs, const Scalar G,
                                               const Mat3 &E, const Mat3 &U1, const Mat3 &U2, const Vec3 &v1,
                                               const Vec3 &v2, const Vec3 &omega1_u, const Vec3 &omega2_u);

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

    for (Index i = 0; i < N; i++) {
        const Vec3 delta_d = dt * v_trans[i] + 0.5 * dt * dt * a_trans[i];
        d_trans[i] += delta_d;
        if (check_energy_balance) {
            delta_d_trans[i] = delta_d;
        }
    }

    // rotations: Simo and Wong algorithm
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
    for (Index i = 0; i < N; i++) {
        const Vec3 a_trans_new = (R_ext_trans[i] - R_int_trans[i]) / M[i];
        v_trans[i] += dt * 0.5 * (a_trans[i] + a_trans_new);
        a_trans[i] = a_trans_new;
    }

    // rotations: Simo and Wong algorithm
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

void calc_forces(const Config &config, const Geometry &geometry, const Borehole &borehole, BeamSystem &beam) {
    const Index N = geometry.get_N();

    zero_internal_and_set_static_forces(N, beam.R_static_trans, beam.R_static_rot, beam.R_int_trans, beam.R_int_rot,
                                        beam.R_ext_trans, beam.R_ext_rot);

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
static void calc_element_inner_forces(const Index ie, const vector<Vec3> &X, const vector<Vec3> &d_trans,
                                      const vector<Quaternion> &d_rot, vector<Vec3> &R_int_trans,
                                      vector<Vec3> &R_int_rot, const Scalar youngs, const Scalar G, const Scalar I_2,
                                      const Scalar I_3, const Scalar A, const Scalar J, const Scalar beta_rayleigh,
                                      const vector<Vec3> &v_trans, const vector<Vec3> &v_rot) {

    // assert(i&e < X.size() - 1);
    const Vec3 &X1 = X[ie];
    const Vec3 &X2 = X[ie + 1];
    const Vec3 &d1 = d_trans[ie];
    const Vec3 &d2 = d_trans[ie + 1];

    DEBUG_ONLY(cout << "X1 " << X1.transpose() << endl; cout << "X2 " << X2.transpose() << endl;);

    DEBUG_ONLY(cout << "d1 " << d1.transpose() << endl; cout << "d2 " << d2.transpose() << endl;);
    const Mat3 T = d_rot[ie].to_matrix();
    const Vec3 &t1 = T.col(0);
    const Vec3 &t2 = T.col(1);
    const Vec3 &t3 = T.col(2);

    const Mat3 U = d_rot[ie + 1].to_matrix();
    const Vec3 &u1 = U.col(0);
    const Vec3 &u2 = U.col(1);
    const Vec3 &u3 = U.col(2);

    DEBUG_ONLY(cout << "U from quat\n"
                    << U << endl;

               cout << "T from quat\n"
                    << T << endl;);

    // calculate the first unit vector
    const Scalar l0 = (X2 - X1).norm();
    const Scalar ln = (X2 + d2 - (X1 + d1)).norm();

    Mat3 E;
    E.col(0) = (X2 + d2 - (X1 + d1)) / ln;
    const Vec3 &e1 = E.col(0);
    assert(is_close(e1.norm(), 1.0));
    DEBUG_ONLY(cout << "e1 = " << e1.transpose() << endl;);
    // Calculate intermediate rotation between triads
    const Mat3 DeltaR = U * T.transpose();
    Quaternion qDeltaR;
    qDeltaR.from_matrix(DeltaR);
    assert(!is_close(qDeltaR.q0, 0.0, 0.1));        // if q0 is 0, gamma_half becomes singular
    const Vec3 gamma_half = qDeltaR.q / qDeltaR.q0; // eq 16.34 divided by 2

    const Mat3 S = skew(gamma_half);
    const Mat3 DeltaR_m = Mat3::Identity() + 1 / (1 + 0.25 * gamma_half.dot(gamma_half)) * (S + 0.5 * S * S);
    assert(is_orthogonal(DeltaR_m));
    const Mat3 R_ = DeltaR_m * T;
    const Vec3 &r1 = R_.col(0);
    const Vec3 &r2 = R_.col(1);
    const Vec3 &r3 = R_.col(2);

    // computing element unit vectors e2 and e3 by rotating the unit vector
    // r1 on to e1.
    E.col(1) = r2 - r2.dot(e1) / (1 + r1.dot(e1)) * (e1 + r1);
    E.col(2) = r3 - r3.dot(e1) / (1 + r1.dot(e1)) * (e1 + r1);
    assert(is_orthogonal(E));
    //    cout << "r1.dot.e1 " << e1.dot(r1) << endl;
    // E.col(1) = r2 - 0.5 * r2.dot(e1) * (e1 + r1);
    // E.col(2) = r3 - 0.5 * r3.dot(e1) * (e1 + r1);
    const Vec3 &e2 = E.col(1);
    const Vec3 &e3 = E.col(2);

    // compute local displacements from global displacements.
    // Taking u_l = ln-l0 is not recommended since subtracting two large
    // numbers may be inaccurate with limited precision. Better to adopt
    // ul = ln - l0 = (ln - l0)*(ln + l0)/(ln + l0) = (ln^2 - l0^2)/(ln + l0)

    const Scalar ul = (ln * ln - l0 * l0) / (ln + l0);
    assert(is_close(ul, ln - l0));

    const Scalar theta_l1 = asin(0.5 * (-t3.dot(e2) + t2.dot(e3)));
    Scalar theta_l2 = asin(0.5 * (-t2.dot(e1) + t1.dot(e2)));
    const Scalar theta_l3 = asin(0.5 * (-t3.dot(e1) + t1.dot(e3)));
    const Scalar theta_l4 = asin(0.5 * (-u3.dot(e2) + u2.dot(e3)));
    const Scalar theta_l5 = asin(0.5 * (-u2.dot(e1) + u1.dot(e2)));
    const Scalar theta_l6 = asin(0.5 * (-u3.dot(e1) + u1.dot(e3)));

    DEBUG_ONLY(cout << "Delta torsion [deg] " << 180 / M_PI * (theta_l4 - theta_l1) << endl;
               cout << "ul " << ul << endl; cout << "theta_l1 " << theta_l1 << endl;
               cout << "theta_l2 " << theta_l2 << endl; cout << "theta_l3 " << theta_l3 << endl;
               cout << "theta_l4 " << theta_l4 << endl; cout << "theta_l5 " << theta_l5 << endl;
               cout << "theta_l6 " << theta_l6 << endl;);

    // cout << "n_glob " << n_glob << endl;
    // assert(theta_l1 == 0);
    // assert(theta_l4 == 0);

#define MAX_ANGLE 90 * M_PI / 180
    assert(abs(theta_l1) < MAX_ANGLE);
    assert(abs(theta_l2) < MAX_ANGLE);
    assert(abs(theta_l3) < MAX_ANGLE);
    assert(abs(theta_l4) < MAX_ANGLE);
    assert(abs(theta_l5) < MAX_ANGLE);
    assert(abs(theta_l6) < MAX_ANGLE);
#undef MAX_ANGLE

    using Mat7_12 = Eigen::Matrix<Scalar, 12, 7>;
    Mat7_12 F_transpose; // Transpose of F where zero rows in F are excluded
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

    const Mat3 A_ = 1.0 / ln * (Mat3::Identity() - e1 * e1.transpose());
    assert(is_close((A_ - A_.transpose()).norm(), 0.0));

    Mat3 L1r2 = 0.5 * r2.dot(e1) * A_ + 0.5 * A_ * r2 * (e1 + r1).transpose();
    Mat3 L1r3 = 0.5 * r3.dot(e1) * A_ + 0.5 * A_ * r3 * (e1 + r1).transpose();
    Mat3 L2r2 = 0.5 * skew(r2) - 0.25 * r2.transpose() * e1 * skew(r1) - 0.25 * skew(r2) * e1 * (e1 + r1).transpose();
    Mat3 L2r3 = 0.5 * skew(r3) - 0.25 * r3.transpose() * e1 * skew(r1) - 0.25 * skew(r3) * e1 * (e1 + r1).transpose();

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

    F_transpose.col(f4) = 1 / (2 * cos(theta_l1)) * (Lr3 * t2 - Lr2 * t3 + h1);
    F_transpose.col(f5) = 1 / (2 * cos(theta_l2)) * (Lr2 * t1 + h2);
    F_transpose.col(f6) = 1 / (2 * cos(theta_l3)) * (Lr3 * t1 + h3);
    F_transpose.col(f10) = 1 / (2 * cos(theta_l4)) * (Lr3 * u2 - Lr2 * u3 + h4);
    F_transpose.col(f11) =
        1 / (2 * cos(theta_l5)) * (Lr2 * u1 + h5); // seems to be an error in the book for f11 and f12. it should be +
                                                   // in front of h6 and h6 (it's + in the paper)
    F_transpose.col(f12) = 1 / (2 * cos(theta_l6)) * (Lr3 * u1 + h6);

    // cout << "F^T\n"
    //      << F.transpose() << endl;

    /*Calculate local internal forces based on linear 3D beam theory*/
    const Vec7 R_int_e_l = calc_element_forces_local(l0, youngs, G, I_2, I_3, A, J, ul, theta_l1, theta_l2, theta_l3,
                                                     theta_l4, theta_l5, theta_l6);

    DEBUG_ONLY(cout << "R_int_e_l:\n" << R_int_e_l << endl;);

    Vec12 R_int_e = F_transpose * R_int_e_l;

    // const Scalar area = M_PI * (ro_e * ro_e - ri_e * ri_e);
    // const Scalar Iy = M_PI / 4 * (ro_e * ro_e * ro_e * ro_e - ri_e * ri_e * ri_e * ri_e);
    // const Scalar Iz = Iy;
    // const Scalar It = 2 * Iy;
    // Vec12 R_int_battini = BattiniBeam::global_internal_forces(ie, X, d_trans, d_rot, area, Iy, Iz, It, youngs, G);
    // assert(R_int_battini.array().isFinite().all());

    // R_int_e = R_int_battini;

    assert(R_int_e.allFinite());
    DEBUG_ONLY(cout << "R_int_e:\n" << R_int_e << endl;);

    DEBUG_ONLY(cout << "R_damp\n" << R_damp << endl;);
    R_int_trans[ie] += R_int_e.segment(0, 3);
    R_int_rot[ie] += R_int_e.segment(3, 3);
    R_int_trans[ie + 1] += R_int_e.segment(6, 3);
    R_int_rot[ie + 1] += R_int_e.segment(9, 3);
    if (beta_rayleigh > 0) {
        const Vec12 R_damp = calc_approx_rayleigh_beta_damping(beta_rayleigh, A, I_2, I_3, J, l0, youngs, G, E, U, T,
                                                               v_trans[ie], v_trans[ie + 1], v_rot[ie], v_rot[ie + 1]);
        R_int_trans[ie] += R_damp.segment(0, 3);
        R_int_rot[ie] += R_damp.segment(3, 3);
        R_int_trans[ie + 1] += R_damp.segment(6, 3);
        R_int_rot[ie + 1] += R_damp.segment(9, 3);
    } else {
        assert(beta_rayleigh == 0.0);
    }

    // DEBUG_ONLY(cout << "R_battini " << R_int_battini << endl;);
    // Vec12 diff = R_int_e - R_int_battini;
    // DEBUG_ONLY(cout << "diff \n"
    //                 << diff << endl;);
    // assert(diff.norm() < 1e-1);
}

static Vec7 calc_element_forces_local(const Scalar l0, const Scalar youngs, const Scalar G, const Scalar I_2,
                                      const Scalar I_3, const Scalar A, const Scalar J, const Scalar ul,
                                      const Scalar theta_1l, const Scalar theta_2l, const Scalar theta_3l,
                                      const Scalar theta_4l, const Scalar theta_5l, const Scalar theta_6l) {

    /*Normal force (F1)
     [[F1], = A*E/l0[[ 1 -1],*[[0],
      [F4]]          [-1  1]]  [ul]]
    */
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
    Scalar M1 = G * J * (theta_1l - theta_4l) / l0;
    Scalar M4 = -M1;

    // M1*=-1;
    // M4*=-1;

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

    Scalar A{}, I_2{}, I_3{}, J{};

#pragma omp parallel for
    for (Index ie = i_first; ie < Ne; ie += 2) {
        geometry.get_cross_section_properties(ie, A, I_2, I_3, J);
        calc_element_inner_forces(ie, X, beam.d_trans, beam.d_rot, beam.R_int_trans, beam.R_int_rot, E, G, I_2, I_3, A,
                                  J, beta_rayleigh, beam.v_trans, beam.v_rot);
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