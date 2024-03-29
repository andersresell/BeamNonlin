
#include "../include/HoleContact.hpp"

void calc_hole_contact_forces(const Config &config, const Index N, const Vec3 *__restrict__ x_hole,
                              Index *__restrict__ hole_index, const Scalar *__restrict__ r_hole,
                              const Scalar *__restrict__ r_outer_string, const Vec3 *__restrict__ X,
                              const Vec3 *__restrict__ d_trans, const Quaternion *__restrict__ d_rot,
                              const Vec3 *__restrict__ v_trans, const Vec3 *__restrict__ v_rot,
                              Vec3 *__restrict__ R_ext_trans, Vec3 *__restrict__ R_ext_rot)
{
    assert(config.contact_enabled);

    update_hole_contact_indices(N, x_hole, hole_index, X, d_trans);

    const Index Ne = N - 1;

    const Scalar mu_s = config.mu_static;
    const Scalar mu_k = config.mu_kinetic;
    const Scalar d_c = config.coloumb_friction_decay;
    const Scalar K_c = config.K_contact;
    const Scalar C_c = config.C_contact;

#pragma omp parallel for
    for (Index i = 0; i < N; i++)
    {
        const Index hi = hole_index[i];
        assert(is_node_within_hole_segment(i, x_hole[hi], x_hole[hi + 1], X[i], d_trans[i]) == 0);
        assert(hi < Ne);
        const Vec3 x = X[i] + d_trans[i];
        const Vec3 &x_hole_A = x_hole[hi];
        const Vec3 &x_hole_B = x_hole[hi + 1];
        const Vec3 t = (x_hole_B - x_hole_A).normalized();
        const Vec3 x_center = x_hole_A + (x - x_hole_A).dot(t) * t;
        const Scalar d_center = (x - x_center).norm();
        const Scalar delta = d_center + r_outer_string[i] - r_hole[i];
        if (delta <= 0.0)
        {
            continue;
        }
        const Scalar dX_segment = dX_element_avg(i, X, N); /*Distance between centroids of two elements*/
        const Vec3 n = (x - x_center).normalized();
        const Vec3 omega = d_rot[i].rotate_vector(v_rot[i]); /*angular velocity in global frame*/

        const Vec3 vc = v_trans[i] + omega.cross(r_outer_string[i] * n); /*Velocity at the contact point*/
        const Scalar vcn = vc.dot(n);
        const Vec3 v_slip = vc - n * vcn;
        const Vec3 t_slip = v_slip.normalized();

        const Scalar mu = mu_k + (mu_s - mu_k) * exp(-d_c * v_slip.norm()); /*Coloumb friction model*/

        const Scalar fn = (K_c * delta + C_c * vcn) * dX_segment; /*Normal force*/
        const Scalar ft = mu * fn;                                /*Tangenital force*/

        const Vec3 fc = -n * fn - t * ft;
        const Vec3 mc = -r_outer_string[i] * n.cross(fc); /*m = r x f*/

        assert(fc.allFinite() && mc.allFinite());

        R_ext_trans[i] += fc;
        R_ext_rot[i] += mc;
    }
}

inline int is_node_within_hole_segment(Index i, const Vec3 &x_hole_A, const Vec3 &x_hole_B,
                                       const Vec3 &X, const Vec3 &d_trans)
{
    const Vec3 t = (x_hole_B - x_hole_A).normalized();
    const Vec3 x = X + d_trans;
    const Scalar l_hole_seg = (x_hole_B - x_hole_A).norm();
    const Scalar dist = (x - x_hole_A).dot(t);
    const bool is_between = dist >= -SMALL_SCALAR && dist <= l_hole_seg + SMALL_SCALAR;
    if (is_between)
    {
        return 0;
    }
    else if (dist < 0)
    {
        return -1;
    }
    else
    {
        assert(dist > 0);
        return 1;
    }
}

inline void update_hole_contact_indices(const Index N, const Vec3 *__restrict__ x_hole,
                                        Index *__restrict__ hole_index, const Vec3 *__restrict__ X,
                                        const Vec3 *__restrict__ d_trans)
{
    const Index Ne = N - 1;
    /*Update hole indices*/
    for (Index i = 0; i < N; i++)
    {
        int hi = hole_index[i];
        int is_between = is_node_within_hole_segment(i, x_hole[hi], x_hole[hi + 1], X[i], d_trans[i]);
        if (is_between == 0)
        {
            continue;
        }
        else if (is_between == -1)
        {
            /*Search backwards*/
            hi--;
            assert(hi >= 0);
            while (is_node_within_hole_segment(i, x_hole[hi], x_hole[hi + 1], X[i], d_trans[i]) != 0)
            {
                assert(is_node_within_hole_segment(i, x_hole[hi], x_hole[hi + 1], X[i], d_trans[i]) == -1);
                hi--;
                assert(hi >= 0);
            }
            hole_index[i] = hi;
        }
        else
        {
            assert(is_between == 1);
            /*Search forwards*/
            hi++;
            assert(hi < Ne);
            while (is_node_within_hole_segment(i, x_hole[hi], x_hole[hi + 1], X[i], d_trans[i]) != 0)
            {
                assert(is_node_within_hole_segment(i, x_hole[hi], x_hole[hi + 1], X[i], d_trans[i]) == 1);
                hi++;
                assert(hi < Ne);
            }
            hole_index[i] = hi;
        }
    }
}

inline Scalar dX_e(Index ie, const Vec3 *X)
{
    return (X[ie + 1] - X[ie]).norm();
}

inline Scalar dX_element_avg(Index i, const Vec3 *X, const Index N)
{
    assert(i < N);
    const Index Ne = N - 1;
    if (i == 0)
        return dX_e(0, X) / 2;
    else if (i == N - 1)
        return dX_e(Ne - 1, X) / 2;
    else
        return (dX_e(i - 1, X) + dX_e(i, X)) / 2;
}