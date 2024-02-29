
#include "../include/Solver.hpp"
#include "../include/InputParser.hpp"

inline void test_quat()
{
    cout << "test quat\n";
    // Scalar thetaz = 15 * M_PI / 180;
    // Scalar thetay = 5 * M_PI / 180;
    // Scalar thetax;

    // Scalar c = cos(thetaz);
    // Scalar s = sin(thetaz);

    // Mat3 Rx, Ry, Rz;
    // Rz << c, -s, 0,
    //     s, c, 0,
    //     0, 0, 1;
    // c = cos(thetay);
    // s = sin(thetay);
    // Ry << c, 0, -s,
    //     0, 1, 0,
    //     s, 0, c;

    // // Ry = Mat3::Identity();

    // Mat3 U = Rz * Ry;

    // assert(is_orthogonal(U));

    // // cout << "U before \n"
    // //      << U << endl;
    // Quaternion q;
    // // cout << "q first \n"
    // //      << q << endl;
    // q.from_matrix(U);

    // // cout << "q second \n"
    // //      << q << endl;
    // // cout << "q norm " << q.norm() << endl;
    // U = q.to_matrix();

    // cout << "q third \n"
    //      << q << endl;

    // cout << "U after \n"
    //      << U << endl;
    Index n_tests = 100;
    for (Index i = 0; i < n_tests; i++)
    {
        Vec3 Theta = Vec3::Random() * 10;

        Vec3 e = Theta.normalized();
        Scalar theta = Theta.norm();
        Mat3 R = Mat3::Identity() + sin(theta) * skew_symmetric(e) + (1 - cos(theta)) * skew_symmetric(e) * skew_symmetric(e);
        assert(is_orthogonal(R));

        Quaternion q;
        q.from_matrix(R);
        Mat3 U = q.to_matrix();

        if ((U - R).norm() > SMALL_SCALAR)
        {
            cout << "Theta " << Theta << endl;
            cout << "Fail\n";
            cout << "U\n"
                 << U << endl;
            cout << "R\n"
                 << R << endl;
        }
    }

    cout << "second test\n";

    for (Index i = 0; i < n_tests; i++)
    {
        Vec3 Theta = Vec3::Random() * 10;
        Vec3 v0 = Vec3::Random() * 100;

        Quaternion q{Theta};
        Vec3 vnq = q.rotate_vector(v0);
        const Mat3 R = q.to_matrix();
        Vec3 vnr = R * v0;
        if ((vnq - vnr).norm() > SMALL_SCALAR)
        {
            cout << "Theta " << Theta << endl;
            cout << "vector rotation failed\n";
            cout << "vnq " << vnq << endl
                 << "vnr " << vnr << endl;
        }

        Vec3 dTheta = Vec3::Random();
        Quaternion dq{dTheta};
        Mat3 dR = dq.to_matrix();

        q.compound_rotate(dTheta);
        Mat3 R_n_q = q.to_matrix();
        Mat3 R_n = dR * R;

        if ((R_n_q - R_n).norm() > SMALL_SCALAR)
        {
            cout << "compoundd rotation failed\n";
        }
    }

    cout << "third test\n";

    for (Index i = 0; i < n_tests; i++)
    {
        Vec3 Theta = Vec3::Random() * 10;
        Vec3 v0 = Vec3::Random() * 100;

        Quaternion q{Theta};
        Vec3 vnq = q.rotate_vector_reversed(v0);
        const Mat3 R = q.to_matrix();
        Vec3 vnr = R.transpose() * v0;
        if ((vnq - vnr).norm() > SMALL_SCALAR)
        {
            cout << "Theta " << Theta << endl;
            cout << "vector reversed rotation failed\n";
            cout << "vnq " << vnq << endl
                 << "vnr " << vnr << endl;
        }
    }

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
