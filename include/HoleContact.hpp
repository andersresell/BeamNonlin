#pragma once
#include "Borehole.hpp"
#include "Config.hpp"
#include "Containers.hpp"
#include "Geometry.hpp"
#include "SolverUtils.hpp"
void calc_hole_contact_forces(const Config &config, const Index N, const Index N_hole, const vector<Vec3> &x_hole,
                              vector<Index> &hole_index, const vector<Scalar> &r_hole,
                              const vector<Scalar> &r_outer_string, const vector<Vec3> &X, const vector<Vec3> &d_trans,
                              const vector<Quaternion> &d_rot, const vector<Vec3> &v_trans, const vector<Vec3> &v_rot,
                              vector<Vec3> &R_ext_trans, vector<Vec3> &R_ext_rot);

template <bool initial_search>
inline void update_hole_contact_indices(const Index N, const Index N_hole, const vector<Vec3> &x_hole,
                                        const vector<Scalar> &r_e_hole, vector<Index> &hole_index,
                                        const vector<Vec3> &X, const vector<Vec3> &d_trans);

inline int node_within_hole_segment(const Index ie_h, const Index N_hole, const vector<Vec3> &x_hole,
                                    const Scalar r_hole, const Vec3 &x);
/*Distance between centroids of two elements*/
inline Scalar dX_element_avg(Index i, const vector<Vec3> &X, const Index N);

void initialize_hole_contact(const Config &config, const Geometry &geometry, const Borehole &borehole,
                             BeamSystem &beam);