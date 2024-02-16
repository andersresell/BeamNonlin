#pragma once
#include "Includes.hpp"
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

enum class BC_Case
{
    NONE,
    CANTILEVER,
    SIMPLY_SUPPORTED
};
