#pragma once
#include "Config.hpp"

class Borehole
{
    vector<Scalar> r_hole_element; /**/
    vector<Vec3> nodes;

public:
    Index get_N_hole_nodes() const { return nodes.size(); }
    Index get_N_hole_elements() const
    {
        assert(nodes.size() == r_hole_element.size() + 1);
        return r_hole_element.size();
    }
    Borehole(const Config &config, Scalar L_string);
};