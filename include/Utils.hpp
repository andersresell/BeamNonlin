#pragma once
#include "Includes.hpp"
#include "Containers.hpp"
#include <chrono>

class Timer
{
    using Clock = std::chrono::high_resolution_clock;
    using Time = std::chrono::_V2::system_clock::time_point;
    Time start_time, end_time;

public:
    void start_counter();
    void stop_counter();
    Scalar get_elapsed_time_sec() const;
    void print_elapsed_time() const;
};

constexpr Scalar STANDARD_GRAVITY = 9.80665;
constexpr Index N_THREADS_MAX = 10;

enum class BC_Case
{
    NONE,
    CANTILEVER,
    SIMPLY_SUPPORTED
};

const map<string, BC_Case> bc_case_from_string{
    {"none", BC_Case::NONE},
    {"cantilever", BC_Case::CANTILEVER},
    {"simply_supported", BC_Case::SIMPLY_SUPPORTED}};

struct PointLoad
{
    Vec3 load_trans;
    Vec3 load_rot;
    Index i;
    PointLoad(const vector<Scalar> R_tmp, Scalar rel_loc, const vector<Vec3> &X);
};

Mat3 triad_from_euler_angles(Scalar theta_x, Scalar theta_y, Scalar theta_z);
