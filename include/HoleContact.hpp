#pragma once
#include "Config.hpp"
#include "Containers.hpp"
// #include "Geometry.hpp"
#include "SolverUtils.hpp"
#include "Borehole.hpp"

void calc_hole_contact_forces(const Config &config, const Index N, const Index Ne_hole,
                              const Vec3 *__restrict__ x_hole,
                              Index *__restrict__ hole_index, const Scalar *__restrict__ r_hole,
                              const Scalar *__restrict__ r_outer_string, const Vec3 *__restrict__ X,
                              const Vec3 *__restrict__ d_trans, const Quaternion *__restrict__ d_rot,
                              const Vec3 *__restrict__ v_trans, const Vec3 *__restrict__ v_rot,
                              Vec3 *__restrict__ R_ext_trans, Vec3 *__restrict__ R_ext_rot);

inline void update_hole_contact_indices(const Index N, const Index Ne_hole, const Vec3 *__restrict__ x_hole,
                                        Index *__restrict__ hole_index, const Vec3 *__restrict__ X,
                                        const Vec3 *__restrict__ d_trans);

inline int node_within_hole_segment(Index i, const Vec3 &x_hole_A, const Vec3 &x_hole_B,
                                    const Vec3 &X, const Vec3 &d_trans);

/*Distance between centroids of two elements*/
inline Scalar dX_element_avg(Index i, const Vec3 *X, const Index N);