#pragma once

#include "Config.hpp"
#include "Containers.hpp"

class Geometry
{
    vector<Vec3> X;
    vector<Scalar> ro;
    vector<Scalar> ri;
    // vector<Scalar> s;

public:
    Geometry(Scalar L0, Scalar N, Scalar D_outer_uniform, Scalar D_inner_uniform);

        const vector<Vec3> &get_X() const { return X; }
    const vector<Scalar> &get_ro() const { return ro; }
    const vector<Scalar> &get_ri() const { return ri; }

    Index get_Ne() const { return X.size() - 1; }
    Index get_N() const { return X.size(); }

    Scalar dx_e(Index ie) const
    {
        assert(ie < get_Ne());
        return (X[ie + 1] - X[ie]).norm();
    }

    Scalar ro_e(Index ie) const
    {
        assert(ie < get_Ne());
        return (ro[ie] + ro[ie + 1]) / 2;
    }
    Scalar ri_e(Index ie) const
    {
        assert(ie < get_Ne());
        return (ri[ie] + ri[ie + 1]) / 2;
    }
    Scalar A_e(Index ie) const
    {
        assert(ie < get_Ne());
        return M_PI * (ro_e(ie) * ro_e(ie) - ri_e(ie) * ri_e(ie));
    }

    Scalar I_e(Index ie) const
    {
        assert(ie < get_Ne());
        Scalar Ro_e_ = ro_e(ie);
        Scalar Ri_e_ = ri_e(ie);
        return M_PI / 4 * (Ro_e_ * Ro_e_ * Ro_e_ * Ro_e_ - Ri_e_ * Ri_e_ * Ri_e_ * Ri_e_);
    }

    // Scalar J_e(Index ie) const { return 2 * I_e(ie); }
};