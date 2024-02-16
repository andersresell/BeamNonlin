#pragma once
#include "Includes.hpp"
#include "Utils.hpp"

struct Config
{
    Config();
    Index n;
    size_t n_max;
    Scalar t_max;
    Scalar dt;
    Scalar t;
    Scalar nu;
    Scalar E;
    Scalar rho;
    Scalar CFL = 0.9;

    Index n_write = 1;
    string base_dir;
    string output_dir;

    BC_Case bc_case = BC_Case::NONE;

    bool save_csv = true;
    bool gravity_enabled = true;

    Scalar get_G() const { return E / (2 * (1 + nu)); }

    size_t get_n_steps() const
    {
        assert(dt > 0.0 && t_max > 0.0);
        return min(n_max, (size_t)(t_max / dt) + 1);
    }
};
