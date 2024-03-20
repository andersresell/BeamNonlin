#pragma once
#include "Config.hpp"
#include "Containers.hpp"
#include "Geometry.hpp"
#include "SolverUtils.hpp"
#include "Borehole.hpp"

namespace BattiniBeam
{
    inline Vec12 r_vector(const Vec3 &x1, const Vec3 &x2);

    inline Mat12_3 G_matrix(const Vec3 &x1, const Vec3 &x2, const Mat3 &U1, const Mat3 &U2);

    inline Mat6_12 P_matrix(const Vec3 &x1, const Vec3 &x2, const Mat3 &U1, const Mat3 &U2);

    inline Mat3 co_rotating_rotation_matrix_current(const Vec3 &x1, const Vec3 &x2, const Mat3 &U1, const Mat3 &U2);

    inline Mat12 E_matrix(const Vec3 &x1, const Vec3 &x2, const Mat3 &U1, const Mat3 &U2);

    inline Vec7 local_deformations(const Vec3 &x1, const Vec3 &x2, const Vec3 &X1, const Vec3 &X2, const Mat3 &U1, const Mat3 &U2);

    inline Mat3 rotation_vector_inverse_T_s(const Index node_id, const Vec3 &x1, const Vec3 &x2, const Vec3 &X1,
                                            const Vec3 &X2, const Mat3 &U1, const Mat3 &U2);

    inline Mat7_12 b_matrix_g(const Vec3 &x1, const Vec3 &x2, const Mat3 &U1, const Mat3 &U2);
    inline Mat7 b_matrix_a(const Vec3 &x1, const Vec3 &x2, const Vec3 &X1,
                           const Vec3 &X2, const Mat3 &U1, const Mat3 &U2);

    inline Mat7 local_stiffness_matrix(const Vec3 &X1, const Vec3 &X2, const Scalar E, const Scalar A,
                                       const Scalar G, const Scalar It,
                                       const Scalar Iy, const Scalar Iz);

    inline Vec7 local_internal_forces(const Vec3 &x1, const Vec3 &x2, const Vec3 &X1,
                                      const Vec3 &X2, const Mat3 &U1, const Mat3 &U2,
                                      const Scalar E, const Scalar A, const Scalar G,
                                      const Scalar It, const Scalar Iy, const Scalar Iz);

    inline Vec7 local_internal_forces_a(const Vec3 &x1, const Vec3 &x2, const Vec3 &X1,
                                        const Vec3 &X2, const Mat3 &U1, const Mat3 &U2,
                                        const Scalar E, const Scalar A, const Scalar G,
                                        const Scalar It, const Scalar Iy, const Scalar Iz);

    Vec12 global_internal_forces(const Index ie, const Vec3 *__restrict__ X, const Vec3 *__restrict__ d_trans,
                                 const Quaternion *__restrict__ d_rot, const Scalar A,
                                 const Scalar Iy, const Scalar Iz, const Scalar It,
                                 const Scalar youngs, const Scalar G);

    void test();
}
