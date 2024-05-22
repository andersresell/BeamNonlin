#pragma once
#include "BattiniBeam.hpp"
#include "Borehole.hpp"
#include "Config.hpp"
#include "Containers.hpp"
#include "Geometry.hpp"
#include "HoleContact.hpp"
#include "SolverUtils.hpp"

void solve(Config &config, Geometry &geometry, const Borehole &borehole);

void set_initial_configuration(const Config &config, vector<Vec3> &X, vector<Vec3> &d_trans, vector<Quaternion> &d_rot);

void calc_dt(Config &config, const Geometry &geometry, const Borehole &borehole);

void check_energy_balance(const Config &config, const BeamSystem &beam_sys);

void calc_forces(const Config &config, const Geometry &geometry, const Borehole &borehole, BeamSystem &beam);

inline void calc_inner_forces(const Config &config, const Geometry &geometry, BeamSystem &beam);

inline void zero_internal_and_set_static_forces(const Index N, const vector<Vec3> &R_static_trans,
                                                const vector<Vec3> &R_static_rot, vector<Vec3> &R_int_trans,
                                                vector<Vec3> &R_int_rot, vector<Vec3> &R_ext_trans,
                                                vector<Vec3> &R_ext_rot);

inline void add_mass_proportional_rayleigh_damping(Index N, Scalar alpha, const vector<Scalar> &M,
                                                   const vector<Vec3> &v_trans, vector<Vec3> &R_int_trans,
                                                   const vector<Vec3> &J_u, const vector<Quaternion> &d_rot,
                                                   const vector<Vec3> &v_rot, vector<Vec3> &R_int_rot);

inline void calc_element_inner_forces(const Index ie, const vector<Vec3> &X, const vector<Vec3> &d_trans,
                                      const vector<Quaternion> &d_rot, vector<Vec3> &R_int_trans,
                                      vector<Vec3> &R_int_rot, const Scalar ri_e, const Scalar ro_e,
                                      const Scalar youngs, const Scalar G, const Scalar beta_rayleigh,
                                      const vector<Vec3> &v_trans, const vector<Vec3> &v_rot);

inline Vec7 calc_element_forces_local(Scalar ri, Scalar ro, Scalar l0, Scalar E, Scalar G, Scalar ul, Scalar theta_1l,
                                      Scalar theta_2l, Scalar theta_3l, Scalar theta_4l, Scalar theta_5l,
                                      Scalar theta_6l);

// inline void step_central_differences(Scalar dt, Index N, Vec3Quat *__restrict__ u, Vec3Vec3 *__restrict__ v,
//                                      const Scalar *__restrict__ M, const Vec3 *__restrict__ J_u,
//                                      const Vec3Vec3 *__restrict__ R_int, const Vec3Vec3 *__restrict__ R_ext,
//                                      bool check_energy_balance, Scalar &W_int, Scalar &W_ext, Scalar &KE);

void calc_static_loads(const Config &config, const Geometry &geometry, vector<Vec3> &R_static_trans,
                       vector<Vec3> &R_static_rot);

inline void set_simple_bc(const Config &config, const Geometry &geometry, BeamSystem &beam);

// // remember to make this inline again
// void velocity_update_partial(Scalar dt, Index N, const Scalar *__restrict__ M, const Vec3 *__restrict__ J_u,
//                              const Vec3 *__restrict__ R_int_trans, const Vec3 *__restrict__ R_int_rot,
//                              const Vec3 *__restrict__ R_ext_trans, const Vec3 *__restrict__ R_ext_rot,
//                              Vec3 *__restrict__ v_trans, Vec3 *__restrict__ v_rot);

// inline void displacement_update(Scalar dt, Index N, Vec3 *__restrict__ v_trans, Vec3 *__restrict__ v_rot,
//                                 Vec3 *__restrict__ d_trans, Quaternion *__restrict__ d_rot);

// inline void calc_delta_d(Scalar dt, Index N, Vec3 *__restrict__ delta_d_trans, Vec3 *__restrict__ delta_d_rot,
//                          const Vec3 *__restrict__ v_trans, const Vec3 *__restrict__ v_rot);
inline void work_update_partial(Index N, const vector<Vec3> &delta_d_trans, const vector<Vec3> &delta_d_rot,
                                const vector<Vec3> &R_int_trans, const vector<Vec3> &R_int_rot,
                                const vector<Vec3> &R_ext_trans, const vector<Vec3> &R_ext_rot, Scalar &W_ext,
                                Scalar &W_int);
inline void kinetic_energy_update(Index N, const vector<Scalar> &M, const vector<Vec3> &J_u,
                                  const vector<Vec3> &v_trans, const vector<Vec3> &v_rot, Scalar &KE);

// inline void rotate_moment_to_body_frame(Index N, const Quaternion *__restrict__ d_rot, Vec3 *__restrict__ R_int_rot,
//                                         Vec3 *__restrict__ R_ext_rot);
void step_explicit_SW(Config &config, const Geometry &geometry, const Borehole &borehole, BeamSystem &beam);

void step_explicit_NMB(Config &config, const Geometry &geometry, const Borehole &borehole, BeamSystem &beam);

void calc_element_forces_local_rotated_TEST(Scalar ri, Scalar ro, Scalar l0, Scalar E, Scalar G, Scalar ul,
                                            Scalar theta_1l, Scalar theta_2l, Scalar theta_3l, Scalar theta_4l,
                                            Scalar theta_5l, Scalar theta_6l, Vec3 &f1, Vec3 &m1, Vec3 &f2, Vec3 &m2);
