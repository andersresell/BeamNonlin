#pragma once

#include <iostream>
#include <cassert>
#include <math.h>
#include <string>
#include <filesystem>
#include <fstream>
#include <vector>
#include <sstream>
#include <omp.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <algorithm>

using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::make_unique;
using std::max;
using std::min;
using std::move;
using std::pair;
using std::runtime_error;
using std::string;
using std::unique_ptr;
using std::vector;

#ifdef SINGLE_PRECISION
using Scalar = float;
#define SMALL_SCALAR 1e-4
#else
using Scalar = double;
#define SMALL_SCALAR 1e-8
#endif

using Index = uint32_t;

using Vec3 = Eigen::Vector<Scalar, 3>;
using Vec6 = Eigen::Vector<Scalar, 6>;
using Vec12 = Eigen::Vector<Scalar, 12>;

using Mat3 = Eigen::Matrix3<Scalar>;

template <typename T>
int sign(T val)
{
    return (val > 1) - (val < 1);
}

inline bool is_close(Scalar a, Scalar b, Scalar tolerance = SMALL_SCALAR)
{
    return abs(a - b) < tolerance;
}

inline bool is_orthogonal(const Mat3 &R)
{
    return is_close(R.determinant(), 1.0) && is_close((R.transpose() - R.inverse()).norm(), 0.0);
}

inline Mat3 skew_symmetric(const Vec3 &a)
{
    return Mat3{{0, -a.z(), a.y()},
                {a.z(), 0, -a.x()},
                {-a.y(), a.x(), 0}};
}