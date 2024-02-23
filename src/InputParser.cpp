#include "../include/InputParser.hpp"

InputParser::InputParser(const string &input_filename) : input_filename{input_filename}
{
    try
    {
        root_node = YAML::LoadFile(input_filename);
    }
    catch (const exception &e)
    {
        throw runtime_error("Error loading input file '" + input_filename + "', " + string(e.what()) +
                            "\nNote that the input file has to be located in the configurations directory.\n");
    }
    cout << "Input file: " << input_filename << endl;
}

void InputParser::parse_config(Config &config) const
{
    try
    {
        parse_yaml_config_options(config);
    }
    catch (exception &e)
    {
        throw runtime_error("Error parsing config options:\n" + string(e.what()));
    }
    cout << "Config options parsed without errors\n";
}

void InputParser::create_geometry(unique_ptr<Geometry> &geometry) const
{
    string root_name = "geometry";
    try
    {
        Scalar L0 = read_required_option<Scalar>(root_name, "L0");
        Scalar N = read_required_option<Scalar>(root_name, "N");
        Scalar D_outer_uniform = read_required_option<Scalar>(root_name, "D_outer_uniform");
        Scalar D_inner_uniform = read_required_option<Scalar>(root_name, "D_inner_uniform");
        geometry = make_unique<Geometry>(L0, N, D_outer_uniform, D_inner_uniform);
    }
    catch (exception &e)
    {
        throw runtime_error("Error parsing and creating geometry object:\n" + string(e.what()));
    }
}

void InputParser::parse_yaml_config_options(Config &config) const
{
    string root_name_setup = "setup";
    string root_name_properties = "properties";

    /*-------------------------------------------------
    ------------------Solver settings------------------
    -------------------------------------------------*/

    config.n_max = read_required_option<size_t>(root_name_setup, "n_max");
    config.t_max = read_required_option<size_t>(root_name_setup, "t_max");
    config.CFL = read_required_option<Scalar>(root_name_setup, "CFL");
    config.save_csv = read_optional_option<bool>(root_name_setup, "save_csv", true);
    config.n_write = read_optional_option<Index>(root_name_setup, "n_write", 1);
    config.bc_case = read_required_enum_option<BC_Case>(root_name_setup, "bc_case", BC_case_from_string);
    config.gravity_enabled = read_required_option<bool>(root_name_setup, "gravity_enabled");
    if (config.gravity_enabled)
    {
        vector<Scalar> gravity_acc;
        /*Specifying user defined gravity or default standard gravity*/
        if (root_node[root_name_setup]["gravity_acc"])
        {
            gravity_acc = read_required_option<vector<Scalar>>(root_name_setup, "gravity_acc");
            if (gravity_acc.size() != 3)
            {
                throw runtime_error("length of gravity_acc option must be 3, but it is " + to_string(gravity_acc.size()) + "\n");
            }
            else
            {
                config.gravity_acc = {gravity_acc[0], gravity_acc[1], gravity_acc[2]};
            }
        }
        else
        {
            config.gravity_acc = {0, 0, -STANDARD_GRAVITY};
        }
    }

    config.n_threads = read_optional_option<Index>(root_name_setup, "n_threads", 1);
    if (config.n_threads > N_THREADS_MAX)
    {
        throw runtime_error(string("number of threads specified is larger than the maximum allowed value.\nn_specified: " +
                                   to_string(config.n_threads) + ", n_max: " + to_string(N_THREADS_MAX)));
    }

    config.nu = read_required_option<Scalar>(root_name_properties, "nu");
    config.E = read_required_option<Scalar>(root_name_properties, "E");
    config.rho = read_required_option<Scalar>(root_name_properties, "rho");
    if (config.CFL >= 1)
    {
        cout << "Warning, CFL >= 1 (CFL = " << config.CFL << ")\n";
    }
    Scalar rho;
    Scalar CFL = 0.9;

    BC_Case bc_case = BC_Case::NONE;

    bool save_csv = true;
    bool gravity_enabled = true;

    Index n_omp_threads = 1;
}

string InputParser::option_not_specified_msg(string root_name, string option_name) const
{
    string msg = "\"" + option_name + "\" not specified in the input file \"" + input_filename + "\"\n";
    msg += "\"with base name " + root_name + "\"";
    return msg;
}