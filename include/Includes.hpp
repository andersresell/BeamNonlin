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
using std::map;
using std::max;
using std::min;
using std::move;
using std::ostream;
using std::pair;
using std::runtime_error;
using std::string;
using std::to_string;
using std::unique_ptr;
using std::vector;

#define PRINT_DEBUG_MESSAGES 1

#ifndef NDEBUG
#ifdef PRINT_DEBUG_MESSAGES
#define DEBUG_ONLY(arg) arg
#endif
#else
#define DEBUG_ONLY(arg)
#endif

#ifndef _OPENMP
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define omp_set_num_threads(num_threads)
#endif

#ifdef SINGLE_PRECISION
using Scalar = float;
#define SMALL_SCALAR 1e-4
#else
using Scalar = double;
// #define SMALL_SCALAR 1e-8
#define SMALL_SCALAR 1e-8
#endif

using Index = uint32_t;

using Vec3 = Eigen::Vector<Scalar, 3>;
using Vec6 = Eigen::Vector<Scalar, 6>;
using Vec7 = Eigen::Vector<Scalar, 7>;
using Vec12 = Eigen::Vector<Scalar, 12>;

using Mat3 = Eigen::Matrix3<Scalar>;

template <typename T>
int sign(T val)
{
    return (val > 1) - (val < 1);
}

inline bool is_close(Scalar a, Scalar b, Scalar tol = SMALL_SCALAR)
{
    const Scalar scaled_tol = tol * max((Scalar)1.0, max(abs(a), abs(b))); // Need to scale the tolerance so that it work for numbers of various sizes
    // assert(abs(b - a) <= scaled_tol);
    return abs(b - a) <= scaled_tol;
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

template <typename T>
inline T sqr(T val)
{
    return val * val;
}

template <typename T>
inline void print_std_vector(const vector<T> &v, string label = "")
{
    if (label != "")
        cout << label << "\n";
    cout << "\n[";
    Index N = v.size();
    for (Index i = 0; i < N; i++)
    {
        cout << " " << v[i];
        if (i < N - 1)
            cout << ",";
        else
            cout << "]\n";
        cout << "\n";
    }
}

template <typename T>
inline void print_raw_array(const T *v, Index N, string label = "")
{
    if (label != "")
        cout << label << "\n";
    cout << "\n[";
    for (Index i = 0; i < N; i++)
    {
        cout << " " << v[i];
        if (i < N - 1)
            cout << ",";
        else
            cout << "]\n";
        cout << "\n";
    }
}
