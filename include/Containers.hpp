#pragma once
#include "Includes.hpp"

struct Quaternion
{
    Scalar q0, q1, q2, q3;

    Mat3 to_matrix() const;
    void from_matrix(const Mat3 &R);
    Scalar norm() const;
    Scalar norm_sqr() const;
};

#include "../src/Containers.inl"