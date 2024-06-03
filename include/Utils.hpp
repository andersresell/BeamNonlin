#pragma once
#include "Containers.hpp"
#include "Includes.hpp"
#include <chrono>

class Timer {
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

enum class BC_Case {
    NONE,
    CANTILEVER,
    SIMPLY_SUPPORTED
};
const map<string, BC_Case> bc_case_from_string{
    {"none", BC_Case::NONE}, {"cantilever", BC_Case::CANTILEVER}, {"simply_supported", BC_Case::SIMPLY_SUPPORTED}};

enum class CrossSectiontype {
    PIPE,
    RECANGLE
};
const map<string, CrossSectiontype> cross_section_type_from_string{{"pipe", CrossSectiontype::PIPE},
                                                                   {"rectangle", CrossSectiontype::RECANGLE}};
enum class CorotationalFormulation {
    CRISFIELD,
    BATTINI
};
const map<string, CorotationalFormulation> corotational_formulation_from_string{
    {"crisfield", CorotationalFormulation::CRISFIELD}, {"battini", CorotationalFormulation::BATTINI}};

struct PointLoad {
    Vec3 load_trans;
    Vec3 load_rot;
    Index i;
    PointLoad(const vector<Scalar> R_tmp, Scalar rel_loc, const vector<Vec3> &X);
};

Mat3 triad_from_euler_angles(Scalar theta_x, Scalar theta_y, Scalar theta_z, string order = "xyz");
