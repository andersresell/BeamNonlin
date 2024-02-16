#pragma once
#include "Containers.hpp"
#include "Geometry.hpp"

struct Vec3Quat
{
    Vec3 trans;
    Quaternion rot;
    Vec3Quat()
    {
        static_assert(sizeof(Vec3Quat) == sizeof(Scalar) * 7);
        trans = {0, 0, 0};
        rot.from_matrix(Mat3::Identity());
    }
};

struct Vec3Vec3
{
    Vec3 trans, rot;
};

struct BeamSystem
{
    vector<Vec3Quat> u; /*Nodal displacements*/
    vector<Vec3Vec3> v;
    vector<Vec3Vec3> R, R_static; /*Nodal forces*/
    vector<Scalar> M_inv;         /*Mass (pre inverted)*/
    vector<Vec3> J_u;             /*Moment of inertia in body frame*/

    BeamSystem(const Config &config, const Geometry &geometry);
};

void save_csv(const Config &config, const Geometry &geometry, const BeamSystem &beam_system);

void create_output_dir(Config &config);