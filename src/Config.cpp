#include "../include/Config.hpp"

Config::Config()
{

    base_dir = std::filesystem::current_path();
    cout << "Simulation directory: " << base_dir << endl;
}