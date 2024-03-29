#pragma once
#include "Config.hpp"

class Borehole
{
    vector<Vec3> X;
    vector<Scalar> r_hole_element; /**/

public:
    const vector<Vec3> &get_X() const { return X; }
    const vector<Scalar> &get_r_hole_element() const { return r_hole_element; }
    Index get_N_hole_nodes() const { return X.size(); }
    Index get_N_hole_elements() const
    {
        assert(X.size() == r_hole_element.size() + 1);
        return r_hole_element.size();
    }
    Borehole(const Config &config, Scalar L_string);
};