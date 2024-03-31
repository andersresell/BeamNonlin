#pragma once
#include "Config.hpp"
#include "Containers.hpp"
#include "Geometry.hpp"
#include "SolverUtils.hpp"
#include "Borehole.hpp"

void calc_hole_contact_forces(const Config &config, const Index N, const Index Ne_hole,
                              const Vec3 *__restrict__ x_hole,
                              Index *__restrict__ hole_index, const Scalar *__restrict__ r_hole,
                              const Scalar *__restrict__ r_outer_string, const Vec3 *__restrict__ X,
                              const Vec3 *__restrict__ d_trans, const Quaternion *__restrict__ d_rot,
                              const Vec3 *__restrict__ v_trans, const Vec3 *__restrict__ v_rot,
                              Vec3 *__restrict__ R_ext_trans, Vec3 *__restrict__ R_ext_rot);

template <bool initial_search>
inline void update_hole_contact_indices(const Index N, const Index N_hole, const Vec3 *__restrict__ x_hole,
                                        const Scalar *__restrict__ r_e_hole,
                                        Index *__restrict__ hole_index, const Vec3 *__restrict__ X,
                                        const Vec3 *__restrict__ d_trans);

inline int node_within_hole_segment(const Index ie_h, const Index N_hole, const Vec3 *__restrict__ x_hole,
                                    const Scalar r_hole, const Vec3 &x);
/*Distance between centroids of two elements*/
inline Scalar dX_element_avg(Index i, const Vec3 *X, const Index N);

void initialize_hole_contact(const Geometry &geometry, const Borehole &borehole, BeamSystem &beam_system);