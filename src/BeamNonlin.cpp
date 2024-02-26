
#include "../include/Solver.hpp"
#include "../include/InputParser.hpp"

inline void test_quat()
{
    cout << "test quat\n";
    Scalar thetaz = 15 * M_PI / 180;
    Scalar thetay = 5 * M_PI / 180;

    Scalar c = cos(thetaz);
    Scalar s = sin(thetaz);

    Mat3 Rz, Ry;
    Rz << c, -s, 0,
        s, c, 0,
        0, 0, 1;
    c = cos(thetay);
    s = sin(thetay);
    Ry << c, 0, -s,
        0, 1, 0,
        s, 0, c;
    Ry = Mat3::Identity();

    Mat3 U = Rz * Ry;

    assert(is_orthogonal(U));

    cout << "U before \n"
         << U << endl;
    Quaternion q;
    cout << "q first \n"
         << q << endl;
    q.from_matrix(U);

    cout << "q second \n"
         << q << endl;
    cout << "q norm " << q.norm() << endl;
    U = q.to_matrix();

    cout << "q third \n"
         << q << endl;

    cout << "U after \n"
         << U << endl;

    exit(0);
}

int main(int argc, char *argv[])
{
    // test_quat();

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
            input_parser.parse_bcs_and_loads(*geometry, config);
        }

        omp_set_num_threads(config.n_threads);
        cout << "Running with " << config.n_threads << " threads\n";
#pragma omp parallel
        {
            printf("Hello from thread %i\n", omp_get_thread_num());
        }

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
