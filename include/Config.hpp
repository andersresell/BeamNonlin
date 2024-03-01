#pragma once
#include "Includes.hpp"
#include "Utils.hpp"

struct Config
{
    Config();
    Index n;
    size_t n_max;
    Scalar t_max;
    Scalar CFL;
    bool save_csv = true;
    Index n_write;
    string base_dir;
    string output_dir;
    BC_Case bc_case = BC_Case::NONE;
    vector<PointLoad> R_point_static;
    bool check_energy_balance;
    Scalar energy_balance_tol;
    bool rayleigh_damping_mass_enabled;
    Scalar alpha_rayleigh;

    bool gravity_enabled = true;
    Vec3 gravity_acc;
    Index n_threads;

    Scalar dt;
    Scalar t;

    Scalar nu;
    Scalar E;
    Scalar rho;

    Scalar get_G() const { return E / (2 * (1 + nu)); }

    size_t get_n_steps() const
    {
        assert(dt > 0.0 && t_max > 0.0);
        return min(n_max, (size_t)(t_max / dt) + 1);
    }
};
