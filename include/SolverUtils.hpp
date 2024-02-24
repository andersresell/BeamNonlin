#pragma once
#include "Containers.hpp"
#include "Geometry.hpp"

struct Vec3Quat
{
    Vec3 trans;
    Quaternion rot;

    Vec3Quat();

    friend ostream &operator<<(ostream &os, const Vec3Quat &rhs);

    static void print_array(vector<Vec3Quat> arr, string label = "",
                            bool print_trans = true, bool print_rot = false);
};

struct Vec3Vec3
{
    Vec3 trans, rot;
    Vec3Vec3() : trans{Vec3::Zero()}, rot{Vec3::Zero()} {}
    void set_zero()
    {
        trans = Vec3::Zero();
        rot = Vec3::Zero();
    }
    friend ostream &operator<<(ostream &os, const Vec3Vec3 &rhs)
    {
        return os << "{" << rhs.trans.transpose() << "}, {" << rhs.rot.transpose() << "}";
    }
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