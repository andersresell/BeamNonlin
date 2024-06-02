#include "../include/InputParser.hpp"

InputParser::InputParser(const string &input_filename) : input_filename{input_filename} {
    try {
        root_node = YAML::LoadFile(input_filename);
    } catch (const exception &e) {
        if (input_filename.size() > 4) {
            string yml_extension = input_filename.substr(input_filename.size() - 4);
            if (yml_extension != ".yml") {
                throw runtime_error("Error loading input file '" + input_filename + "'\n" +
                                    "Don't forget .yml extension in name");
            }
        }

        throw runtime_error("Error loading input file '" + input_filename + "', " + string(e.what()) +
                            "\nNote that the input file has to be located in the configurations directory.\n");
    }
    cout << "Input file: " << input_filename << endl;
}

void InputParser::parse_config(Config &config) const {
    try {
        parse_yaml_config_options(config);
    } catch (exception &e) {
        throw runtime_error("Error parsing config options:\n" + string(e.what()));
    }
    cout << "Config options parsed without errors\n";
}

void InputParser::create_geometry(unique_ptr<Geometry> &geometry) const {

    Scalar D_outer_uniform;
    Scalar D_inner_uniform;
    Scalar h2_uniform;
    Scalar h3_uniform;

    string root_name = "geometry";
    try {
        Scalar L0 = read_required_option<Scalar>(root_name, "L0");
        Scalar N = read_required_option<Scalar>(root_name, "N");

        CrossSectiontype cross_section_type = read_optional_enum_option<CrossSectiontype>(
            root_name, "cross_section_type", cross_section_type_from_string, CrossSectiontype::PIPE);

        if (cross_section_type == CrossSectiontype::PIPE) {
            D_outer_uniform =
                read_required_option<Scalar>(root_name, "D_outer_uniform", "Must be specified for pipe cross section");
            D_inner_uniform =
                read_required_option<Scalar>(root_name, "D_inner_uniform", "Must be specified for pipe cross section");
        } else {
            assert(cross_section_type == CrossSectiontype::RECANGLE);
            h2_uniform =
                read_required_option<Scalar>(root_name, "h2_uniform", "Must be specified for rectangle cross section");
            h3_uniform =
                read_required_option<Scalar>(root_name, "h3_uniform", "Must be specified for rectangle cross section");
        }
        geometry =
            make_unique<Geometry>(L0, N, D_outer_uniform, D_inner_uniform, h2_uniform, h3_uniform, cross_section_type);
    } catch (exception &e) {
        throw runtime_error("Error parsing and creating geometry object:\n" + string(e.what()));
    }
}

void InputParser::parse_yaml_config_options(Config &config) const {
    string root_name_setup = "setup";
    string root_name_properties = "properties";

    config.n_max = read_required_option<size_t>(root_name_setup, "n_max");
    config.t_max = read_required_option<Scalar>(root_name_setup, "t_max");
    config.CFL = read_required_option<Scalar>(root_name_setup, "CFL");
    config.save_csv = read_optional_option<bool>(root_name_setup, "save_csv", true);
    // config.user_defined_force = read_optional_option<bool>(root_name_setup, "user_defined_force", false);
    // printf("Warning: Using user defined force, remember to include the source file when building, or else nothing
    // will "
    //        "happen\n");
    config.n_write = read_optional_option<Index>(root_name_setup, "n_write", 1);
    config.check_energy_balance = read_optional_option<bool>(root_name_setup, "check_energy_balance", false);
    if (config.check_energy_balance) {
        config.energy_balance_tol = read_optional_option<Scalar>(root_name_setup, "energy_balance_tol", 0.01);
    }
    config.rayleigh_damping_enabled = read_optional_option<bool>(root_name_setup, "rayleigh_damping_enabled", false);
    if (config.rayleigh_damping_enabled) {
        config.alpha_rayleigh = read_required_option<Scalar>(root_name_setup, "alpha_rayleigh");
        config.beta_rayleigh = read_required_option<Scalar>(root_name_setup, "beta_rayleigh");
    } else {
        config.alpha_rayleigh = 0;
        config.beta_rayleigh = 0;
    }
    config.borehole_included = read_optional_option<bool>(root_name_setup, "borehole_included", false);

    config.corotational_beam_formulation = read_optional_enum_option<CorotationalBeamFormulation>(
        root_name_setup, "corotational_beam_formulation", corotational_beam_formulation_from_string,
        CorotationalBeamFormulation::CRISFIELD);

    config.gravity_enabled = read_required_option<bool>(root_name_setup, "gravity_enabled");
    if (config.gravity_enabled) {
        vector<Scalar> gravity_acc;
        /*Specifying user defined gravity or default standard gravity*/
        if (root_node[root_name_setup]["gravity_acc"]) {
            gravity_acc = read_required_option<vector<Scalar>>(root_name_setup, "gravity_acc");
            if (gravity_acc.size() != 3) {
                throw runtime_error("length of gravity_acc option must be 3, but it is " +
                                    to_string(gravity_acc.size()) + "\n");
            } else {
                config.gravity_acc = {gravity_acc[0], gravity_acc[1], gravity_acc[2]};
            }
        } else {
            config.gravity_acc = {0, 0, -STANDARD_GRAVITY};
        }
    }
    config.contact_enabled = read_optional_option<bool>(root_name_setup, "contact_enabled", false);
    if (config.contact_enabled) {
        config.K_contact = read_required_option<Scalar>(root_name_properties, "K_contact");
        config.C_contact = read_required_option<Scalar>(root_name_properties, "C_contact");
        config.coloumb_friction_decay = read_required_option<Scalar>(root_name_properties, "coloumb_friction_decay");
        config.mu_static = read_required_option<Scalar>(root_name_properties, "mu_static");
        config.mu_kinetic = read_required_option<Scalar>(root_name_properties, "mu_kinetic");
    }

    config.n_threads = read_optional_option<Index>(root_name_setup, "n_threads", 1);
    if (config.n_threads > N_THREADS_MAX) {
        throw runtime_error(
            string("number of threads specified is larger than the maximum allowed value.\nn_specified: " +
                   to_string(config.n_threads) + ", n_max: " + to_string(N_THREADS_MAX)));
    }

    config.nu = read_required_option<Scalar>(root_name_properties, "nu");
    config.E = read_required_option<Scalar>(root_name_properties, "E");
    config.rho = read_required_option<Scalar>(root_name_properties, "rho");
    if (config.CFL >= 1) {
        cout << "Warning, CFL >= 1 (CFL = " << config.CFL << ")\n";
    }
}
void InputParser::parse_bcs_and_loads(const Geometry &geometry, Config &config) const {
    string root_name_bc = "bc";
    string root_name_loads = "loads";
    try {
        config.bc_case = read_required_enum_option<BC_Case>(root_name_bc, "case", bc_case_from_string);
        if (config.bc_case == BC_Case::CANTILEVER) {

            if (root_node[root_name_bc]["orientation_base_euler_angles_xyz_deg"]) {
                vector<Scalar> orient_euler_deg =
                    read_required_option<vector<Scalar>>(root_name_bc, "orientation_base_euler_angles_xyz_deg");
                if (orient_euler_deg.size() != 3) {
                    throw runtime_error("Wrong size specified for the option: "
                                        "\"orientation_base_euler_angles_xyz_deg\", the size must be 3");
                }
                const Scalar theta_x = orient_euler_deg[0] * M_PI / 180;
                const Scalar theta_y = orient_euler_deg[1] * M_PI / 180;
                const Scalar theta_z = orient_euler_deg[2] * M_PI / 180;
                config.bc_orientation_base.from_matrix(triad_from_euler_angles(theta_x, theta_y, theta_z));
            } else {
                config.bc_orientation_base.from_matrix(Mat3::Identity());
            }
        } else {
            assert(false);
        }
    } catch (exception &e) {
        throw runtime_error("Error parsing boundary conditions:\n" + string(e.what()));
    }

    try {
        if (root_node[root_name_loads]["point_loads"]) {
            YAML::Node point_loads = root_node[root_name_loads]["point_loads"];
            vector<Scalar> R_point_tmp;
            Scalar rel_loc;
            for (const auto &point_load : point_loads) {
                config.point_loads_rel_to_base_orientation =
                    read_optional_option<bool>(root_name_loads, "point_loads_rel_to_base_orientation", false);
                if (point_load["R"] && point_load["rel_loc"]) {
                    R_point_tmp = point_load["R"].as<vector<Scalar>>();
                    if (R_point_tmp.size() != 6) {
                        throw runtime_error("Wrong size for point load vector R, must be 6, but is " +
                                            to_string(R_point_tmp.size()));
                    }
                    rel_loc = point_load["rel_loc"].as<Scalar>();
                    if (rel_loc < 0 || rel_loc > 1) {
                        throw runtime_error("0 <= rel_loc <= 1 must be fulfilled, but rel_loc = " + to_string(rel_loc));
                    }
                    config.R_point_static.push_back(PointLoad{R_point_tmp, rel_loc, geometry.get_X()});
                } else {
                    throw runtime_error("Point load incorrectly specified\n");
                }
            }
        }
    } catch (exception &e) {
        throw runtime_error("Error parsing loads:\n" + string(e.what()));
    }
}

string InputParser::option_not_specified_msg(string root_name, string option_name, string extra_msg) const {
    string msg = "\"" + option_name + "\" not specified in the input file \"" + input_filename + "\"\n";
    msg += "\"with base name \"" + root_name + "\"";
    if (!extra_msg.empty())
        msg += "\n" + extra_msg;
    return msg;
}
