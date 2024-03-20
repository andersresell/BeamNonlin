#pragma once
#include "Config.hpp"
#include "Containers.hpp"
#include "Geometry.hpp"
#include "SolverUtils.hpp"
#include "Borehole.hpp"

namespace BattiniBeam
{
    inline Vec12 r_vector(const Vec3 &x1, const Vec3 &x2)
    {
        Vec3 r_1 = x2 - x1;
        r_1 /= r_1.norm();

        Vec12 r_vec = Vec12::Zero();
        r_vec.segment(0, 3) = -r_1;
        r_vec.segment(6, 3) = r_1;
        return r_vec;
    }

    inline Mat12_3 G_matrix(const Vec3 &x1, const Vec3 &x2, const Mat3 &U1, const Mat3 &U2)
    {
        const Mat3 reference_rotation = Mat3::Identity();
        const Mat3 current_rotation = co_rotating_rotation_matrix_current(x1, x2, U1, U2);

        const Vec3 q_1 = U1 * reference_rotation * Vec3{0, 1, 0};
        const Vec3 q_2 = U2 * reference_rotation * Vec3{0, 1, 0};

        const Vec3 q_mean = 0.5 * (q_1 + q_2);

        const Vec3 q_i = ((current_rotation.transpose()) * q_mean);
        const Vec3 q_1i = ((current_rotation.transpose()) * q_1);
        const Vec3 q_2i = ((current_rotation.transpose()) * q_1);

        const Scalar eta = q_i[0] / q_i[1];
        const Scalar eta_11 = q_1i[0] / q_i[1];
        const Scalar eta_12 = q_1i[1] / q_i[1];
        const Scalar eta_21 = q_2i[0] / q_i[1];
        const Scalar eta_22 = q_2i[1] / q_i[1];

        const Scalar l_n = (x2 - x1).norm();
        Mat3_12 G_transposed = Mat3_12::Zero();
        G_transposed(0, 2) = eta / l_n;
        G_transposed(0, 3) = eta_12 / 2.0;
        G_transposed(0, 4) = -eta_11 / 2.0;
        G_transposed(0, 8) = -eta / l_n;
        G_transposed(0, 9) = eta_22 / 2.0;
        G_transposed(0, 10) = -eta_21 / 2.0;
        return G_transposed.transpose();
    }

    inline Mat6_12 P_matrix(const Vec3 &x1, const Vec3 &x2, const Mat3 &U1, const Mat3 &U2)
    {
        Mat6_12 P = Mat6_12::Zero();
        P.block<3, 3>(0, 0) = Mat3::Identity();
        P.block<3, 3>(3, 3) = Mat3::Identity();

        Mat12_3 G_transposed = G_matrix(x1, x2, U1, U2).transpose();
        Mat6_12 p_temp = Mat6_12::Zero();
        p_temp.block<12, 3>(0, 0) = G_transposed;
        p_temp.block<6, 12>(3, 0) = G_transposed;

        P = P - p_temp;
        return P;
    }

    inline Mat3 co_rotating_rotation_matrix_current(const Vec3 &x1, const Vec3 &x2, const Mat3 &U1, const Mat3 &U2)
    {
        Mat3 R = Mat3::Zero();

        const Mat3 reference_rotation = Mat3::Identity();

        const Vec3 r_1 = (x2 - x1).normalized();

        const Vec3 q_1 = U1 * reference_rotation * Vec3{0, 1, 0};
        const Vec3 q_2 = U2 * reference_rotation * Vec3{0, 1, 0};

        const Vec3 q_mean = 0.5 * (q_1 + q_2);

        const Vec3 r_3 = r_1.cross(q_mean).normalized();
        const Vec3 r_2 = r_3.cross(r_1).normalized();

        R.col(0) = r_1;
        R.col(1) = r_2;
        R.col(2) = r_3;
        return R;
    }

    inline Mat12 E_matrix(const Vec3 &x1, const Vec3 &x2, const Mat3 &U1, const Mat3 &U2)
    {
        Mat3 R_r = co_rotating_rotation_matrix_current(x1, x2, U1, U2);
        Mat12 E = Mat12::Zero();
        E.block<3, 3>(0, 0) = R_r;
        E.block<3, 3>(3, 3) = R_r;
        E.block<3, 3>(6, 6) = R_r;
        E.block<3, 3>(9, 9) = R_r;
        return E;
    }

    inline Vec7 local_deformations(const Vec3 &x1, const Vec3 &x2, const Vec3 &X1, const Vec3 &X2, const Mat3 &U1, const Mat3 &U2)
    {
        Vec7 disp_l = Vec7::Zero();
        disp_l[0] = ((x2 - X2) - (x1 - X1)).norm();

        const Mat3 R_co_rotating_reference = Mat3::Identity();
        const Mat3 R_co_rotating_current = co_rotating_rotation_matrix_current(x1, x2, U1, U2);

        const Mat3 R_local_node_1 = R_co_rotating_current.transpose() * U1 * R_co_rotating_reference;
        const Mat3 R_local_node_2 = R_co_rotating_current.transpose() * U2 * R_co_rotating_reference;

        Scalar local_rot_norm_node_1 = acos(R_local_node_1.trace() - 1.0) / 2.0;
        Scalar local_rot_norm_node_2 = acos(R_local_node_2.trace() - 1.0) / 2.0;

        if (local_rot_norm_node_1 != 0.00)
        {
            Mat3 log_R_local_node_1 = (local_rot_norm_node_1 / (2.0 * sin(local_rot_norm_node_1))) * (R_local_node_1 - (R_local_node_1.transpose()));
            disp_l[1] = log_R_local_node_1(2, 1);
            disp_l[2] = log_R_local_node_1(0, 2);
            disp_l[3] = log_R_local_node_1(1, 0);
        }
        if (local_rot_norm_node_2 != 0.00)
        {
            Mat3 log_R_local_node_2 = (local_rot_norm_node_2 / (2.0 * sin(local_rot_norm_node_2))) * (R_local_node_2 - (R_local_node_2.transpose()));
            disp_l[4] = log_R_local_node_2(2, 1);
            disp_l[5] = log_R_local_node_2(0, 2);
            disp_l[6] = log_R_local_node_2(1, 0);
        }

        return disp_l;
    }

    inline Mat3 rotation_vector_inverse_T_s(const Index node_id, const Vec3 &x1, const Vec3 &x2, const Vec3 &X1,
                                            const Vec3 &X2, const Mat3 &U1, const Mat3 &U2)
    {
        Vec3 disp_l;
        if (node_id == 1)
        {
            disp_l = local_deformations(x1, x2, X1, X2, U1, U2).segment(1, 3);
        }
        else if (node_id == 2)
        {
            disp_l = local_deformations(x1, x2, X1, X2, U1, U2).segment(4, 3);
        }
        else
        {
            printf("node nr: %i not known", node_id);
            exit(1);
        }

        const Scalar phi_norm = disp_l.norm();
        Mat3 T_inv = Mat3::Identity();
        if (phi_norm > 0.0)
        {
            const Vec3 u_vec = disp_l / phi_norm;
            const Mat3 phi_tilde{{0.0, -disp_l[2], disp_l[1]},
                                 {disp_l[2], 0.0, -disp_l[0]},
                                 {-disp_l[1], disp_l[0], 0.0}};

            T_inv = ((phi_norm / 2.0) / (tan(phi_norm / 2.0))) * Mat3::Identity();
            T_inv += (1.0 - ((phi_norm / 2.0) / (tan(phi_norm / 2.0)))) * u_vec * u_vec.transpose();
            T_inv -= 0.5 * phi_tilde;
        }
        assert(is_orthogonal(T_inv));
        return T_inv;
    }

    inline Mat7_12 b_matrix_g(const Vec3 &x1, const Vec3 &x2, const Mat3 &U1, const Mat3 &U2)
    {
        const Vec12 r_vec = r_vector(x1, x2);
        const Mat6_12 P_E_t = P_matrix(x1, x2, U1, U2) * E_matrix(x1, x2, U1, U2).transpose();
        Mat7_12 b_mat_g = Mat7_12::Zero();
        b_mat_g.row(0) = r_vec;
        b_mat_g.block<6, 12>(1, 0) = P_E_t;
        return b_mat_g;
    }

    inline Mat7 b_matrix_a(const Vec3 &x1, const Vec3 &x2, const Vec3 &X1,
                           const Vec3 &X2, const Mat3 &U1, const Mat3 &U2)
    {
        Mat7 b = Mat7::Zero();
        b(0, 0) = 1.0;
        b.block<3, 3>(1, 1) = rotation_vector_inverse_T_s(1, x1, x2, X1, X2, U1, U2);
        b.block<3, 3>(4, 4) = rotation_vector_inverse_T_s(2, x1, x2, X1, X2, U1, U2);
        return b;
    }

    inline Mat7 local_stiffness_matrix(const Vec3 &X1, const Vec3 &X2, const Scalar E, const Scalar A,
                                       const Scalar G, const Scalar It,
                                       const Scalar Iy, const Scalar Iz)
    {
        Mat7 C = Mat7::Zero();
        C(0, 0) = E * A;
        C(1, 1) = G * It;
        C(4, 4) = G * It;
        C(1, 4) = -G * It;
        C(4, 1) = -G * It;
        C(2, 2) = 4.0 * E * Iz;
        C(2, 5) = 2.0 * E * Iz;
        C(3, 3) = 4.0 * E * Iy;
        C(3, 6) = 2.0 * E * Iy;
        C(5, 5) = 4.0 * E * Iz;
        C(5, 2) = 2.0 * E * Iz;
        C(6, 6) = 4.0 * E * Iy;
        C(6, 3) = 2.0 * E * Iy;
        C /= (X2 - X1).norm();
        return C;
    }

    inline Vec7 local_internal_forces(const Vec3 &x1, const Vec3 &x2, const Vec3 &X1,
                                      const Vec3 &X2, const Mat3 &U1, const Mat3 &U2,
                                      const Scalar E, const Scalar A, const Scalar G,
                                      const Scalar It, const Scalar Iy, const Scalar Iz)
    {
        const Vec7 loc_def = local_deformations(x1, x2, X1, X2, U1, U2);
        const Mat7 K_l = local_stiffness_matrix(X1, X2, E, A, G, It, Iy, Iz);
        return K_l * loc_def;
    }

    inline Mat7 local_internal_forces_a(const Vec3 &x1, const Vec3 &x2, const Vec3 &X1,
                                        const Vec3 &X2, const Mat3 &U1, const Mat3 &U2,
                                        const Scalar E, const Scalar A, const Scalar G,
                                        const Scalar It, const Scalar Iy, const Scalar Iz)
    {
        return b_matrix_a(x1, x2, X1, X2, U1, U2).transpose() *
               local_internal_forces(x1, x2, X1, X2, U1, U2, E, A, G, It, Iy, Iz);
    }

    inline void calc_element_inner_forces(const Index ie, const Vec3 *__restrict__ X, const Vec3 *__restrict__ d_trans,
                                          const Quaternion *__restrict__ d_rot, Vec3 *__restrict__ R_int_trans,
                                          Vec3 *__restrict__ R_int_rot, const Scalar ri, const Scalar ro, const Scalar youngs,
                                          const Scalar G, const Scalar beta_rayleigh, const Vec3 *__restrict__ v_trans,
                                          const Vec3 *__restrict__ v_rot)
    {
        const Vec3 X1 = X[ie];
        const Vec3 X2 = X[ie + 1];
        const Vec3 x1 = X1 + d_trans[ie];
        const Vec3 x2 = X2 + d_trans[ie + 1];
        const Mat3 U1 = d_rot[ie].to_matrix();
        const Mat3 U2 = d_rot[ie + 1].to_matrix();

        const Scalar A = M_PI * (ro * ro - ri * ri);

        const Scalar Iy = M_PI / 4 * (ro * ro * ro * ro - ri * ri * ri * ri);
        const Scalar Iz = Iy;
        const Scalar It = 2 * Iy;

        const Vec12 f_int_g = b_matrix_g(x1, x2, U1, U2).transpose() * local_internal_forces_a(x1, x2, X1, X2, U1, U2, youngs, A, G, It, Iy, Iz);
    }
}
