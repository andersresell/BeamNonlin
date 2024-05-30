#pragma once
#include "Borehole.hpp"
#include "Config.hpp"
#include "Containers.hpp"
#include "Geometry.hpp"
#include "HoleContact.hpp"
#include "SolverUtils.hpp"

void step_explicit(Config &config, const Geometry &geometry, const Borehole &borehole, BeamSystem &beam);
