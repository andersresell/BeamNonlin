#pragma once
#include "Borehole.hpp"
#include "Containers.hpp"
#include "Geometry.hpp"

extern Index n_glob;
extern Scalar t_glob;
// struct BeamSystem
// {
//     vector<Vec3Quat> u; /*Nodal displacements*/
//     vector<Vec3Vec3> v;
//     vector<Vec3Vec3> R_int, R_ext, R_static; /*Nodal forces*/
//     vector<Scalar> M;                        /*Lumped mass */
//     vector<Vec3> J_u;                        /*Moment of inertia in body frame*/
//     Scalar W_int, W_ext, KE;                 /*Variables for energy balance check*/
//     vector<Vec3Vec3> delta_u;                /*Used if checking energy balance*/

//     BeamSystem(const Config &config, const Geometry &geometry);
// };
struct BeamSystem {
    vector<Vec3> d_trans;
    vector<Quaternion> d_rot;
    vector<Vec3> v_trans;
    vector<Vec3> a_trans;
    vector<Vec3> v_rot;                                    // angular velocity in body frame
    vector<Vec3> a_rot;                                    // angular acceleration in body frame
    vector<Vec3> L_rot;                                    // Angular impulse
    vector<Vec3> m_rot;                                    // Moment at t_{n+1/2}
    vector<Vec3> R_int_trans, R_ext_trans, R_static_trans; /*Nodal translational forces*/

    vector<Vec3> R_int_rot, R_ext_rot, R_static_rot; /*Nodal rotational forces (moments)*/
    vector<Scalar> M;                                /*Lumped mass */
    vector<Vec3> J_u;                                /*Moment of inertia in body frame*/
    Scalar W_int, W_ext, KE;                         /*Variables for energy balance check*/
    vector<Vec3> delta_d_trans;                      /*Used if checking energy balance*/
    vector<Vec3> delta_d_rot;
    vector<Index> hole_index; /*Used to refer to the index of the nearest hole if contact is enabled*/

    BeamSystem(const Config &config, const Geometry &geometry);
};

void save_csv(const Config &config, const Geometry &geometry, const BeamSystem &beam);

void write_borehole(const Config &config, const Borehole &borehole);

void create_output_dir(Config &config);
