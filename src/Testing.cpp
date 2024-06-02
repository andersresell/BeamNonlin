#include "../include/Testing.hpp"
#include "../include/BattiniBeam.hpp"
#include "../include/CrisfieldBeam.hpp"
#include "../include/Numerics.hpp"

void test_inner_forces() {
    printf("Testing inner forces\n");
    // massive pipe
    Scalar L0 = 10;
    Scalar Diameter = 0.1;
    Scalar youngs = 210000000000.0;
    Scalar nu = 0.3;
    Scalar A = M_PI / 4 * powi<2>(Diameter);
    Scalar I_2 = M_PI / 64 * powi<4>(Diameter);
    Scalar I_3 = I_2;
    Scalar J = 2 * I_2;
    Scalar G = youngs / (2.0 * (1 + nu));

    vector<Vec3> X;
    X.push_back(Vec3{0, 0, 0});
    X.push_back(Vec3{L0, 0, 0});

    printf("Case: clamped in both ends + vertical disp\n");
    // Clamped in both end and vertical disp
    Scalar Delta_z = -0.01;
    vector<Vec3> d_trans;
    d_trans.push_back(Vec3{0, 0, 0});
    d_trans.push_back(Vec3{0, 0, Delta_z});

    vector<Quaternion> d_rot;
    Quaternion q;
    q.from_matrix(Mat3::Identity());
    d_rot.push_back(q);
    d_rot.push_back(q);

    vector<Vec3> v_trans;
    vector<Vec3> v_rot;

    Vec12 R_int_crisfield =
        CrisfieldBeam::calc_element_inner_forces(0, X, d_trans, d_rot, youngs, G, I_2, I_3, A, J, 0, v_trans, v_rot);

    Vec12 R_int_battini = BattiniBeam::global_internal_forces(0, X, d_trans, d_rot, A, I_2, I_3, J, youngs, G);

    cout << "Crisfield:\n"
         << R_int_crisfield << endl
         << "\nBattini:\n"
         << R_int_battini << endl
         << "\ndiff:\n"
         << (R_int_crisfield - R_int_battini) << endl;

    assert(is_close(1.0 / powi<3>(2), 1.0 / 8));
    Mat4 kEB_z = Mat4{{12, 6 * L0, -12, 6 * L0},
                      {6 * L0, 4 * L0 * L0, -6 * L0, 2 * L0 * L0},
                      {-12, -6 * L0, 12, -6 * L0},
                      {6 * L0, 2 * L0 * L0, -6 * L0, 4 * L0 * L0}} *
                 youngs * I_2 / powi<3>(L0);
    Vec4 loc_disp_z = Vec4{0, 0, Delta_z, 0};
    Vec4 R_loc = kEB_z * loc_disp_z;

    cout << "N battini " << R_int_battini[0] << ", N crisfield " << R_int_crisfield[0] << endl;

    cout << "V1 theory " << R_loc[0] << ", V1 Battini " << R_int_battini[2] << ", V1 Crisfield " << R_int_crisfield[2]
         << endl
         << "M1 theory " << R_loc[1] << ", M1 Battini " << R_int_battini[4] << ", M1 Crisfield " << R_int_crisfield[4]
         << endl
         << "V2 theory " << R_loc[2] << ", V2 Battini " << R_int_battini[8] << ", V2 Crisfield " << R_int_crisfield[8]
         << endl
         << "M2 theory " << R_loc[3] << ", M2 Battini " << R_int_battini[10] << ", M2 Crisfield " << R_int_crisfield[10]
         << endl;

    exit(0);
}

void test_quat() {
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
    for (Index i = 0; i < n_tests; i++) {
        Vec3 Theta = Vec3::Random() * 10;

        Vec3 e = Theta.normalized();
        Scalar theta = Theta.norm();
        Mat3 R = Mat3::Identity() + sin(theta) * skew(e) + (1 - cos(theta)) * skew(e) * skew(e);
        assert(is_orthogonal(R));

        Quaternion q;
        q.from_matrix(R);
        Mat3 U = q.to_matrix();

        if ((U - R).norm() > SMALL_SCALAR) {
            cout << "Theta " << Theta << endl;
            cout << "Fail\n";
            cout << "U\n" << U << endl;
            cout << "R\n" << R << endl;
        }
    }

    cout << "second test\n";

    for (Index i = 0; i < n_tests; i++) {
        Vec3 Theta = Vec3::Random() * 10;
        Vec3 v0 = Vec3::Random() * 100;

        Quaternion q{Theta};
        Vec3 vnq = q.rotate_vector(v0);
        const Mat3 R = q.to_matrix();
        Vec3 vnr = R * v0;
        if ((vnq - vnr).norm() > SMALL_SCALAR) {
            cout << "Theta " << Theta << endl;
            cout << "vector rotation failed\n";
            cout << "vnq " << vnq << endl << "vnr " << vnr << endl;
        }

        Vec3 dTheta_u = Vec3::Random();
        Quaternion dq{dTheta_u};
        Mat3 dR_u = dq.to_matrix();

        // q.compound_rotate(dTheta);
        q.exponential_map_body_frame(dTheta_u);

        Mat3 R_n_q = q.to_matrix();
        Mat3 R_n = R * dR_u;

        if ((R_n_q - R_n).norm() > SMALL_SCALAR) {
            cout << "compoundd rotation failed\n";
        }
    }

    cout << "third test\n";

    for (Index i = 0; i < n_tests; i++) {
        Vec3 Theta = Vec3::Random() * 10;
        Vec3 v0 = Vec3::Random() * 100;

        Quaternion q{Theta};
        Vec3 vnq = q.rotate_vector_reversed(v0);
        const Mat3 R = q.to_matrix();
        Vec3 vnr = R.transpose() * v0;
        if ((vnq - vnr).norm() > SMALL_SCALAR) {
            cout << "Theta " << Theta << endl;
            cout << "vector reversed rotation failed\n";
            cout << "vnq " << vnq << endl << "vnr " << vnr << endl;
        }
    }

    exit(0);
}