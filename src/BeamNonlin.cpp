
#include "../include/Solver.hpp"
#include "../include/InputParser.hpp"

int main(int argc, char *argv[])
{
    try
    {
        cout << "\n\n////////////////////////////////////////////\n"
             << "////////// Corotational Beam Solver ///////////\n"
             << "///////////////////////////////////////////////\n\n";
        if (argc != 2)
        {
            throw runtime_error("Specify the yaml input file\n");
        }

        string input_file = argv[1];

        // Scalar L0 = 2;
        // Scalar N = 2;
        // Scalar D_outer = 0.1;
        // Scalar D_inner = 0.05;

        Config config{};
        // config.n_max = 5;
        // config.n_write = 500 * 8;
        // config.save_csv = true;
        // config.t_max = 1000000;
        // config.nu = 0.3;
        // config.E = 200 * pow(10, 9);
        // config.rho = 7000;
        // config.CFL = 0.9;
        // config.bc_case = BC_Case::CANTILEVER;
        // config.n_omp_threads = 4;

        // L0 = 1000;
        // N = 300;
        // config.CFL = 0.1;
        // config.n_max = 50000 * 8;
        // config.save_csv = false;
        unique_ptr<Geometry> geometry;
        {
            InputParser input_parser{input_file};
            input_parser.parse_config(config);
            input_parser.create_geometry(geometry);
        }

        omp_set_num_threads(config.n_threads);
        cout << "Running with " << config.n_threads << " threads\n";
        printf("Hello from thread %i\n", omp_get_thread_num());

        solve(config, *geometry);
    }
    catch (exception &e)
    {
        cerr << "Exception caught:\n"
             << e.what() << endl;
        exit(EXIT_FAILURE);
    }
}
