#pragma once
#include "Borehole.hpp"
#include "Config.hpp"
#include "Containers.hpp"
#include "Geometry.hpp"
#include "SolverUtils.hpp"

namespace BattiniBeam {
Vec12 r_vector(const Vec3 &x1, const Vec3 &x2);

Mat12_3 G_matrix(const Vec3 &x1, const Vec3 &x2, const Mat3 &U1, const Mat3 &U2);

Mat6_12 P_matrix(const Vec3 &x1, const Vec3 &x2, const Mat3 &U1, const Mat3 &U2);

Mat3 co_rotating_rotation_matrix_current(const Vec3 &x1, const Vec3 &x2, const Mat3 &U1, const Mat3 &U2);

Mat12 E_matrix(const Vec3 &x1, const Vec3 &x2, const Mat3 &U1, const Mat3 &U2);

Vec7 local_deformations(const Vec3 &x1, const Vec3 &x2, const Vec3 &X1, const Vec3 &X2, const Mat3 &U1, const Mat3 &U2);

Mat3 rotation_vector_inverse_T_s(const Index node_id, const Vec3 &x1, const Vec3 &x2, const Vec3 &X1, const Vec3 &X2,
                                 const Mat3 &U1, const Mat3 &U2);

Mat7_12 b_matrix_g(const Vec3 &x1, const Vec3 &x2, const Mat3 &U1, const Mat3 &U2);
Mat7 b_matrix_a(const Vec3 &x1, const Vec3 &x2, const Vec3 &X1, const Vec3 &X2, const Mat3 &U1, const Mat3 &U2);

Mat7 local_stiffness_matrix(const Vec3 &X1, const Vec3 &X2, const Scalar E, const Scalar A, const Scalar G,
                            const Scalar It, const Scalar Iy, const Scalar Iz);

Vec7 local_internal_forces(const Vec3 &x1, const Vec3 &x2, const Vec3 &X1, const Vec3 &X2, const Mat3 &U1,
                           const Mat3 &U2, const Scalar E, const Scalar A, const Scalar G, const Scalar It,
                           const Scalar Iy, const Scalar Iz);

Vec7 local_internal_forces_a(const Vec3 &x1, const Vec3 &x2, const Vec3 &X1, const Vec3 &X2, const Mat3 &U1,
                             const Mat3 &U2, const Scalar E, const Scalar A, const Scalar G, const Scalar It,
                             const Scalar Iy, const Scalar Iz);

Vec12 global_internal_forces(const Index ie, const vector<Vec3> &X, const vector<Vec3> &d_trans,
                             const vector<Quaternion> &d_rot, const Scalar A, const Scalar Iy, const Scalar Iz,
                             const Scalar It, const Scalar youngs, const Scalar G);

void test();
} // namespace BattiniBeam
