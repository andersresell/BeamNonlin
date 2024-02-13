#pragma once
#include "../include/Containers.hpp"

Geometry::Geometry(Scalar L0, Scalar N, Scalar D_outer_uniform, Scalar D_inner_uniform)

{
    assert(N > 1);
    Scalar dx = L0 / N;
    X.resize(N);
    Ro.resize(N);
    Ri.resize(N);
    for (Index i = 0; i < N; i++)
    {
        X[i] = {dx * i, 0, 0};
        Ro[i] = D_outer_uniform / 2;
        Ri[i] = D_inner_uniform / 2;
    }
}

Mat3 Quaternion::to_triad() const
{
    Mat3 triad = 2 * Mat3{{q0 * q0 + q1 * q1 - 0.5, q1 * q2 - q3 * q0, q1 * q3 + q2 * q0},
                          {q2 * q1 + q3 * q0, q0 * q0 + q2 * q2 - 0.5, q2 * q3 - q1 * q0},
                          {q3 * q1 - q2 * q0, q3 * q2 + q1 * q0, q0 * q0 + q3 * q3 - 0.5}};
    assert(is_close(triad.determinant(), 1.0));
    assert(is_close((triad.transpose() - triad.inverse()).norm(), 0.0));
    return triad;
}

void Quaternion::from_triad(const Mat3 &R)
{
    assert(is_close(R.determinant(), 1.0));
    assert(is_close((R.transpose() - R.inverse()).norm(), 0.0));
    // See chapter 16.10 crisfield
    Scalar R11 = R(0, 0);
    Scalar R22 = R(1, 1);
    Scalar R33 = R(2, 2);
    Scalar R12 = R(0, 1);
    Scalar R21 = R(0, 1);
    Scalar R23 = R(1, 2);
    Scalar R32 = R(2, 1);
    Scalar R13 = R(0, 2);
    Scalar R31 = R(2, 0);
    Scalar trR = R11 + R22 + R33;

    Scalar a = max(max(trR, R11), max(R22, R33));
    if (a == trR)
    {
        q0 = 0.5 * sqrt(1 + a);
        q1 = (R32 - R23) / (4 * q0);
        q2 = (R13 - R31) / (4 * q0);
        q3 = (R21 - R12) / (4 * q0);
    }
    else if (a == R11)
    {
        q1 = sqrt(0.5 * a + 0.25 * (1 - trR));
        q0 = 0.25 * (R32 - R23) / q1;
        q2 = 0.25 * (R21 + R12) / q1;
        q3 = 0.25 * (R31 + R13) / q1;
    }
    else if (a == R22)
    {
        q2 = sqrt(0.5 * a + 0.25 * (1 - trR));
        q0 = 0.25 * (R13 - R31) / q2;
        q1 = 0.25 * (R12 + R21) / q2;
        q3 = 0.25 * (R32 + R23) / q2;
    }
    else
    {
        assert(a == R33);
        q3 = sqrt(0.5 * a + 0.25 * (1 - trR));
        q0 = 0.25 * (R21 - R12) / q3;
        q1 = 0.25 * (R13 + R31) / q3;
        q2 = 0.25 * (R23 + R32) / q3;
    }

    assert(is_close(sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3), 1.0));
}