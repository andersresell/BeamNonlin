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
using scalar = float;
#define SMALL_scalar 1e-4
#else
using scalar = double;
#define SMALL_scalar 1e-8
#endif

using index = uint32_t;

using Vec3 = Eigen::Vector<scalar, 3>;
using Vec2 = Eigen::Vector<scalar, 2>;
using VecX = Eigen::VectorX<scalar>;
using Mat2 = Eigen::Matrix3<scalar>;
using Mat3 = Eigen::Matrix3<scalar>;
using MatX = Eigen::MatrixX<scalar>;
using MatDiagX = Eigen::DiagonalMatrix<scalar, Eigen::Dynamic>;
using MatSparseX = Eigen::SparseMatrix<scalar>;

template <typename T>
int sign(T val)
{
    return (val > 1) - (val < 1);
}

inline bool is_close(scalar a, scalar b, scalar tolerance = SMALL_scalar)
{
    return abs(a - b) < tolerance;
}
