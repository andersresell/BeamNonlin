#include "../include/Utils.hpp"

Scalar Timer::get_elapsed_time_sec() const
{
    using std::chrono::duration_cast;
    using Microseconds = std::chrono::microseconds;
    size_t total_microseconds = duration_cast<Microseconds>(end_time - start_time).count();
    return total_microseconds / 1000000.0;
}

void Timer::print_elapsed_time() const
{
    cout << "Elapsed time = " << get_elapsed_time_sec() << " sec\n";
};

void Timer::start_counter() { start_time = Clock::now(); }
void Timer::stop_counter() { end_time = Clock::now(); }

PointLoad::PointLoad(const vector<Scalar> R_tmp, Scalar rel_loc, const vector<Vec3> &X)
{
    assert(R_tmp.size() == 6);

    this->load_trans = {R_tmp[0], R_tmp[1], R_tmp[2]};
    this->load_rot = {R_tmp[3], R_tmp[4], R_tmp[5]};

    /*Finding the nearest node index corresponding to the relative location*/

    const Index N = X.size();
    assert(N >= 2);
#ifndef NDEBUG
    for (const auto &Xi : X)
    {
        assert(Xi.y() == 0 && Xi.z() == 0); // only straigt beams aligned along x axis allowed for now
    }
#endif
    Scalar L0 = (X[X.size() - 1] - X[0]).norm();
    Scalar x = rel_loc * L0;
    assert(x >= 0 && x <= L0);

    for (Index i = 0; i < N; i++)
    {
        if (x <= X[i].x())
        {
            this->i = i;
            break;
        }
        else if (i == N - 1)
        {
            assert(false);
        }
    }
}

Mat3 triad_from_euler_angles(Scalar theta_x, Scalar theta_y, Scalar theta_z)
{
    Scalar c, s;
    c = cos(theta_x);
    s = sin(theta_x);
    Mat3 Ux = Mat3{{1, 0, 0},
                   {0, c, -s},
                   {0, s, c}};
    c = cos(theta_y);
    s = sin(theta_y);
    Mat3 Uy = Mat3{{c, 0, s},
                   {0, 1, 0},
                   {-s, 0, c}};
    c = cos(theta_z);
    s = sin(theta_z);
    Mat3 Uz = Mat3{{c, -s, 0},
                   {s, c, 0},
                   {0, 0, 1}};

    Mat3 U = Uz * Uy * Ux;
    assert(is_orthogonal(U));
    return U;
}