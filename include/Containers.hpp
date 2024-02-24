#pragma once
#include "Includes.hpp"

struct Quaternion
{
    Scalar q0, q1, q2, q3;

    Mat3 to_matrix() const;
    void from_matrix(const Mat3 &R);
    Scalar norm() const;
    Scalar norm_sqr() const;
    void normalize();
    friend ostream &operator<<(ostream &os, const Quaternion &rhs)
    {
        os << "q0=" << rhs.q0 << "\nq1=" << rhs.q1 << "\nq2=" << rhs.q2 << "\nq3=" << rhs.q3 << endl;
    }
};

#include "../src/Containers.inl"