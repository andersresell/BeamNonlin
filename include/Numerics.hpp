#pragma once
#include "Borehole.hpp"
#include "Config.hpp"
#include "Geometry.hpp"
#include "SolverUtils.hpp"

void step_explicit(Config &config, const Geometry &geometry, const Borehole &borehole, BeamSystem &beam);

void calc_initial_accelerations(const Config &config, const Geometry &geometry, const Borehole &borehole,
                                BeamSystem &beam);