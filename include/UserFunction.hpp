#pragma once
#include "SolverUtils.hpp"



/*Declaration for user defined force. Implement this and to run a specific simulation*/
void add_user_defined_external_forces(const Config &config, const Geometry &geometry, const Borehole &borehole,
                                      const vector<Vec3> &d_trans, const vector<Quaternion> &d_rot,
                                      vector<Vec3> &R_ext_trans, vector<Vec3> &R_ext_rot);

void print_user_defined_external_forces_message();