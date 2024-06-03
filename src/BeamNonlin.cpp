
#include "../include/Borehole.hpp"
#include "../include/Includes.hpp"
#include "../include/InputParser.hpp"
#include "../include/Solver.hpp"
#include "../include/Testing.hpp"

int main(int argc, char *argv[]) {
    // test_quat();
    // test_inner_forces();
    // test_quaternion_performance_comparison(100000000, Vec3{1, 1, 1}, Vec3{1, 0, 0});

    try {

        cout << "\n\n"
             << "///////////////////////////////////////////////\n"
             << "////////// Corotational Beam Solver ///////////\n"
             << "///////////////////////////////////////////////\n\n";
        if (argc != 2) {
            throw runtime_error("Specify the yaml input file\n");
        }

        string input_file = argv[1];
        Config config{};

        unique_ptr<Geometry> geometry;
        {
            InputParser input_parser{input_file};
            input_parser.parse_config(config);
            input_parser.create_geometry(geometry);
            input_parser.parse_bcs_and_loads(*geometry, config);
        }
        unique_ptr<Borehole> borehole = make_unique<Borehole>(config, geometry->get_L0());

        omp_set_num_threads(config.n_threads);
        cout << "Running with " << config.n_threads << " threads\n";
#pragma omp parallel
        { printf("Hello from thread %i\n", omp_get_thread_num()); }
        solve(config, *geometry, *borehole);
        exit(0);
    } catch (exception &e) {
        cerr << "Exception caught:\n" << e.what() << endl;
        exit(EXIT_FAILURE);
    }
}
