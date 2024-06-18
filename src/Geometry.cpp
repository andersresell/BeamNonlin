#include "../include/Geometry.hpp"

Geometry::Geometry(Scalar L0, Index N, Scalar D_outer_uniform, Scalar D_inner_uniform, Scalar h2_uniform,
                   Scalar h3_uniform, CrossSectiontype cross_section_type)
    : cross_section_type{cross_section_type} {
    assert(N > 1);
    Scalar dx = L0 / (N - 1);
    X.resize(N);
    for (Index i = 0; i < N; i++) {
        X[i] = {dx * i, 0, 0};
    }
    if (cross_section_type == CrossSectiontype::PIPE) {
        ro.resize(N, D_outer_uniform / 2);
        ri.resize(N, D_inner_uniform / 2);
    } else {
        assert(cross_section_type == CrossSectiontype::RECANGLE);
        h2.resize(N, h2_uniform);
        h3.resize(N, h3_uniform);
    }
}
