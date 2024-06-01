#include "../include/UserFunction.hpp"

// void print_user_defined_external_forces_message() {
// #ifdef USER_DEFINED_EXTERNAL_FORCE_FILE
//     printf("--------------------------------------------------------\n"
//            "WARNING! USER DEFINED EXTERNAL FORCES ACTIVATED.\n"
//            "--------------------------------------------------------\n");
// #endif
// }

void add_user_defined_external_forces(const Config &config, const Geometry &geometry, const Borehole &borehole,
                                      const vector<Vec3> &d_trans, const vector<Quaternion> &d_rot,
                                      vector<Vec3> &R_ext_trans, vector<Vec3> &R_ext_rot) {
    // dummy function
}