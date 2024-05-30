#pragma once

#include "Config.hpp"
#include "Containers.hpp"

class Geometry {

    vector<Vec3> X;
    vector<Scalar> ro; /*Nodal cross section outer radius*/
    vector<Scalar> ri; /*Nodal cross section inner radius*/
    vector<Scalar> h2; /*Nodal cross section x2-height*/
    vector<Scalar> h3; /*Nodal cross section x3-height*/
    CrossSectiontype cross_section_type;

  public:
    Geometry(Scalar L0, Scalar N, Scalar D_outer_uniform, Scalar D_inner_uniform, Scalar h2_uniform, Scalar h3_uniform,
             CrossSectiontype cross_section_type);

    const vector<Vec3> &get_X() const { return X; }
    vector<Vec3> &get_X() { return X; }
    const vector<Scalar> &get_ro() const { return ro; }
    const vector<Scalar> &get_ri() const { return ri; }

    Index get_Ne() const { return X.size() - 1; }
    Index get_N() const { return X.size(); }
    Scalar get_L0() const {
        assert(X[0].isZero() && X[X.size() - 1].y() == 0 && X[X.size() - 1].z() == 0);
        return (X[X.size() - 1] - X[0]).norm();
    }

    Scalar dx_e(Index ie) const {
        assert(ie < get_Ne());
        return (X[ie + 1] - X[ie]).norm();
    }

    Scalar ro_e(Index ie) const {
        assert(cross_section_type == CrossSectiontype::PIPE);
        assert(ie < get_Ne());
        return (ro[ie] + ro[ie + 1]) / 2;
    }
    Scalar ri_e(Index ie) const {
        assert(cross_section_type == CrossSectiontype::PIPE);
        assert(ie < get_Ne());
        return (ri[ie] + ri[ie + 1]) / 2;
    }

    Scalar h2_e(Index ie) const {
        assert(cross_section_type == CrossSectiontype::RECANGLE);
        assert(ie < get_Ne());
        return (h2[ie] + h2[ie + 1]) / 2;
    }
    Scalar h3_e(Index ie) const {
        assert(cross_section_type == CrossSectiontype::RECANGLE);
        assert(ie < get_Ne());
        return (h3[ie] + h3[ie + 1]) / 2;
    }

    Scalar A_e(const Index ie) const {
        assert(ie < get_Ne());
        if (cross_section_type == CrossSectiontype::PIPE) {
            return M_PI * (powi<2>(ro_e(ie)) - powi<2>(ri_e(ie)));
        } else {
            assert(cross_section_type == CrossSectiontype::RECANGLE);
            return h2_e(ie) * h3_e(ie);
        }
    }

    void get_cross_section_properties(const Index ie, Scalar &A, Scalar &I_2, Scalar &I_3, Scalar &J) const {
        assert(ie < get_Ne());
        A = A_e(ie);
        if (cross_section_type == CrossSectiontype::PIPE) {
            const Scalar ro_ = ro_e(ie);
            const Scalar ri_ = ri_e(ie);
            I_2 = M_PI / 4 * (powi<4>(ro_) - powi<4>(ri_));
            I_3 = I_2;
            J = 2 * I_2;
        } else {
            assert(cross_section_type == CrossSectiontype::RECANGLE);
            const Scalar h2_ = h2_e(ie);
            const Scalar h3_ = h3_e(ie);
            I_2 = h2_ * powi<3>(h3_) / 12;
            I_2 = h3_ * powi<3>(h2_) / 12;
            J = h2_ * h3_ * (h2_ * h2_ + h3_ * h3_) / 12;
        }
    }
};
