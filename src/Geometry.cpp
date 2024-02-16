#include "../include/Geometry.hpp"

Geometry::Geometry(Scalar L0, Scalar N, Scalar D_outer_uniform, Scalar D_inner_uniform)

{
    assert(N > 1);
    Scalar dx = L0 / (N - 1);
    X.resize(N);
    ro.resize(N);
    ri.resize(N);
    for (Index i = 0; i < N; i++)
    {
        X[i] = {dx * i, 0, 0};
        ro[i] = D_outer_uniform / 2;
        ri[i] = D_inner_uniform / 2;
    }
}
