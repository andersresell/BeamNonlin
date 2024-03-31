#pragma once
#include "Config.hpp"

class Borehole
{
    vector<Vec3> x;                /*Hole nodes*/
    vector<Scalar> r_hole_element; /*Radii of hole elements*/

public:
    const vector<Vec3> &get_x() const { return x; }
    const vector<Scalar> &get_r_hole_element() const { return r_hole_element; }
    Index get_N_hole_nodes() const { return x.size(); }
    Index get_N_hole_elements() const
    {
        assert(x.size() == r_hole_element.size() + 1);
        return r_hole_element.size();
    }
    Borehole(const Config &config, Scalar L_string);
    Scalar dx(Index ie) const
    {
        assert(ie < x.size() - 1);
        return (x[ie + 1] - x[ie]).norm();
    }
};