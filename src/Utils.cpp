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

    this->load.trans = {R_tmp[0], R_tmp[1], R_tmp[2]};
    this->load.rot = {R_tmp[3], R_tmp[4], R_tmp[5]};

    /*Finding the nearest node index corresponding to the relative location*/

    const Index N = X.size();
    assert(N >= 2);
    for (const auto &Xi : X)
    {
        assert(Xi.y() == 0 && Xi.z() == 0); // only straigt beams aligned along x axis allowed for now
    }
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