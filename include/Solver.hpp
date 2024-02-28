#pragma once
#include "Config.hpp"
#include "Containers.hpp"
#include "Geometry.hpp"
#include "SolverUtils.hpp"

void solve(Config &config, const Geometry &geometry);

void calc_dt(Config &config, const Geometry &geometry);

static Index n_glob;

inline void assemble(const Config &config, const Geometry &geometry, BeamSystem &beam_sys);

inline void calc_element_inner_forces(Index ie, const Vec3 *__restrict__ X, const Vec3Quat *__restrict__ d,
                                      Vec3Vec3 *__restrict__ R, Scalar ri_e, Scalar ro_e, Scalar E, Scalar G);

inline Vec7 calc_nodal_forces_local(Scalar ri, Scalar ro, Scalar l0, Scalar E, Scalar G, Scalar ul,
                                    Scalar theta_1l, Scalar theta_2l, Scalar theta_3l, Scalar theta_4l,
                                    Scalar theta_5l, Scalar theta_6l);

inline void step_central_differences(Scalar dt,  Vec3Quat*__restrict__ u, Vec3Vec3*__restrict__ v, const Scalar*__restrict__ M_inv,
                                     const Vec3 *__restrict__ J_u, Vec3Vec3*__restrict__ R, const Vec3Vec3*__restrict__ R_static);

inline void calc_static_loads(const Config &config, const Geometry &geometry, vector<Vec3Vec3> &R_static);

inline void set_simple_bc(BC_Case bc_case, const Geometry &geometry, BeamSystem &beam_system);

#include "../src/Solver.inl"
