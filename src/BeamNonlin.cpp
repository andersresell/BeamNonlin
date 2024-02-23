
#include "../include/Solver.hpp"
#include "../include/InputParser.hpp"

int main(int argc, char *argv[])
{
    try
    {
        cout << "\n\n"
             << "///////////////////////////////////////////////\n"
             << "////////// Corotational Beam Solver ///////////\n"
             << "///////////////////////////////////////////////\n\n";
        if (argc != 2)
        {
            throw runtime_error("Specify the yaml input file\n");
        }

        string input_file = argv[1];

        Config config{};

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
        exit(0);
    }
    catch (exception &e)
    {
        cerr << "Exception caught:\n"
             << e.what() << endl;
        exit(EXIT_FAILURE);
    }
}
