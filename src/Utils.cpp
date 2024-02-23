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