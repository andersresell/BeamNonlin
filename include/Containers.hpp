#pragma once
#include "Includes.hpp"

struct Quaternion
{
    Scalar q0, q1, q2, q3;

    Mat3 to_triad() const;

    void from_triad(const Mat3 &R);
};

class Geometry
{
    vector<Vec3> X;
    vector<Scalar> Ro;
    vector<Scalar> Ri;
    // vector<Scalar> s;

public:
    Geometry(Scalar L0, Scalar N, Scalar D_outer_uniform, Scalar D_inner_uniform);

    Index get_Ne() const { return X.size() - 1; }
    Index get_N() const { return X.size(); }

    Scalar dx_e(Index ie) const
    {
        assert(ie < get_Ne());
        return (X[ie + 1] - X[ie]).norm();
    }

    Scalar Ro_e(Index ie) const
    {
        assert(ie < get_Ne());
        return (Ro[ie] * Ro[ie + 1]) / 2;
    }
    Scalar Ri_e(Index ie) const
    {
        assert(ie < get_Ne());
        return (Ri[ie] * Ri[ie + 1]) / 2;
    }
    Scalar A_e(Index ie)
    {
        assert(ie < get_Ne());
        return M_PI * (Ro_e(ie) * Ro_e(ie) - Ri_e(ie) * Ri_e(ie));
    }
};

#include "../src/Containers.inl"