#include "../include/Borehole.hpp"

Borehole::Borehole(const Config &config, Scalar L_string)
{
    if (!config.borehole_included)
    {
        return;
    }
    // Just hardcoding in an L ish shape for now
    Scalar r_uniform = 0.2;
    Index N_hole = 20;
    nodes.resize(N_hole);
    r_hole_element.resize(N_hole - 1, r_uniform);
    Scalar L_hole = L_string * 1.3;
    Scalar dx = L_hole / N_hole;
    for (Index i = 0; i < N_hole; i++)
    {
        Scalar x = dx * i;
        nodes[i] = {x, cos(x / L_hole * M_PI), cos(2 * x / L_hole * M_PI)};
    }
}