#pragma once
#include "Config.hpp"
#include "Containers.hpp"
#include "SolverRuntime.hpp"

void solve(Config &config, const Geometry &geometry);

void calc_dt(Config &config, const Geometry &geometry);