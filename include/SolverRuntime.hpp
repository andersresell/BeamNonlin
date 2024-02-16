#pragma once
#include "Utils.hpp"
#include "SolverUtils.hpp"
#include "Geometry.hpp"

inline void assemble(const Config &config, const Geometry &geometry, BeamSystem &beam_sys);

inline void calc_element_contribution(Index ie, const vector<Vec3> &X, vector<Vec3Quat> &u,
                                      vector<Vec3Vec3> &R, Scalar ri, Scalar ro, Scalar E, Scalar G);

inline Vec12 calc_nodal_forces_local(Scalar ri, Scalar ro, Scalar l0, Scalar E, Scalar G, Scalar ul,
                                     Scalar theta_l1, Scalar theta_l2, Scalar theta_l3, Scalar theta_l4,
                                     Scalar theta_l5, Scalar theta_l6);

inline void step_central_differences(Scalar dt, vector<Vec3Quat> &u, vector<Vec3Vec3> &v, vector<Scalar> M_inv,
                                     const vector<Vec3> &J_u, const vector<Vec3Vec3> &R, const vector<Vec3Vec3> &R_static);

inline void calc_static_loads(const Config &config, const Geometry &geometry, vector<Vec3Vec3> &R_static);

inline void set_simple_bc(BC_Case bc_case, const Geometry &geometry, BeamSystem &beam_system);

#include "../src/SolverRuntime.inl"
