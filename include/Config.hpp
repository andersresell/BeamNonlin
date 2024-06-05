#pragma once
#include "Includes.hpp"
#include "Utils.hpp"

struct Config {
    Config();
    Index n;
    size_t n_max;
    Scalar t_max;
    Scalar CFL;
    // bool user_defined_force;
    Index n_write;
    string base_dir;
    string output_dir;
    BC_Case bc_case = BC_Case::NONE;
    Quaternion bc_orientation_base;
    vector<PointLoad> R_point_static;
    Scalar energy_balance_tol;
    Scalar alpha_rayleigh;
    Scalar beta_rayleigh;
    CorotationalFormulation corotational_formulation;

    Vec3 gravity_acc;
    Index n_threads;

    Scalar dt;
    Scalar t;

    Scalar nu;
    Scalar E;
    Scalar rho;

    Scalar mu_static;
    Scalar mu_kinetic;
    Scalar coloumb_friction_decay;
    Scalar K_contact;
    Scalar C_contact;

    bool save_csv = true;
    bool point_loads_rel_to_base_orientation;
    bool check_energy_balance;
    bool rayleigh_damping_enabled;
    bool borehole_included;
    bool gravity_enabled = true;
    bool contact_enabled;

    Scalar get_G() const { return E / (2 * (1 + nu)); }

    size_t get_n_steps() const {
        assert(dt > 0.0 && t_max > 0.0);
        return min(n_max, (size_t)(t_max / dt) + 1);
    }
};
