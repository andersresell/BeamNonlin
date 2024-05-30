
#include "../include/HoleContact.hpp"

void calc_hole_contact_forces(const Config &config, const Index N, const Index N_hole, const vector<Vec3> &x_hole,
                              vector<Index> &hole_index, const vector<Scalar> &r_hole,
                              const vector<Scalar> &r_outer_string, const vector<Vec3> &X, const vector<Vec3> &d_trans,
                              const vector<Quaternion> &d_rot, const vector<Vec3> &v_trans, const vector<Vec3> &v_rot,
                              vector<Vec3> &R_ext_trans, vector<Vec3> &R_ext_rot) {
    assert(config.contact_enabled);

    update_hole_contact_indices<false>(N, N_hole, x_hole, r_hole, hole_index, X, d_trans);

    const Scalar mu_s = config.mu_static;
    const Scalar mu_k = config.mu_kinetic;
    const Scalar d_c = config.coloumb_friction_decay;
    const Scalar K_c = config.K_contact;
    const Scalar C_c = config.C_contact;

#pragma omp parallel for
    for (Index i = 0; i < N; i++) {
        const Index hi = hole_index[i];
        assert(hi < N_hole - 1);
        const Vec3 x = X[i] + d_trans[i];
        assert(node_within_hole_segment(hi, N_hole, x_hole, r_hole[hi], x) == 0);
        const Vec3 &x_hole_A = x_hole[hi];
        const Vec3 &x_hole_B = x_hole[hi + 1];
        const Vec3 t = (x_hole_B - x_hole_A).normalized();
        const Vec3 x_center = x_hole_A + (x - x_hole_A).dot(t) * t;
        const Scalar d_center = (x - x_center).norm();
        const Scalar delta = d_center + r_outer_string[i] - r_hole[hi];
        if (delta <= 0.0) {
            continue;
        }

        const Scalar dX_segment = dX_element_avg(i, X, N); /*Distance between centroids of two elements*/
        const Vec3 n = (x - x_center).normalized();
        const Vec3 omega = d_rot[i].rotate_vector(v_rot[i]);             /*angular velocity in global frame*/
        const Vec3 vc = v_trans[i] + omega.cross(r_outer_string[i] * n); /*Velocity at the contact point*/
        const Scalar vcn = vc.dot(n);
        const Vec3 v_slip = vc - n * vcn;
        const Vec3 t_slip = v_slip.normalized();

        const Scalar mu = mu_k + (mu_s - mu_k) * exp(-d_c * v_slip.norm()); /*Coloumb friction model*/

        const Scalar fn = (K_c * delta + C_c * vcn) * dX_segment; /*Normal force*/
        const Scalar ft = mu * fn;                                /*Tangenital force*/
        const Vec3 fc = -n * fn - t_slip * ft;
        const Vec3 mc = r_outer_string[i] * n.cross(fc); /*m = r x f*/
        assert(fc.allFinite() && mc.allFinite());

        R_ext_trans[i] += fc;
        R_ext_rot[i] += mc;
    }
}

inline int node_within_hole_segment(const Index ie_h, const Index N_hole, const vector<Vec3> &x_hole,
                                    const Scalar r_hole, const Vec3 &x) {
    const Index Ne_hole = N_hole - 1;
    assert(ie_h < Ne_hole);

    const Vec3 &x_hole_A = x_hole[ie_h];
    const Vec3 &x_hole_B = x_hole[ie_h + 1];

    const Vec3 t = (x_hole_B - x_hole_A).normalized();
    const Scalar l_hole_seg = (x_hole_B - x_hole_A).norm();
    const Scalar dist = (x - x_hole_A).dot(t);

    const bool is_between = dist >= -SMALL_SCALAR && dist <= l_hole_seg + SMALL_SCALAR;

#define MAX_ALLOWABLE_HOLE_PAIR_ANGLE 10 * M_PI / 180
    if (is_between) {
        return 0;
    } else if (dist < 0) {
        if (ie_h == 0) {
            return -1;
        } else {
            assert(ie_h > 0);
            const Vec3 &x_hole_prev = x_hole[ie_h - 1];
            const Vec3 t_prev = (x_hole_A - x_hole_prev).normalized();
            const Scalar theta = acos(t.dot(t_prev)) + M_PI;
            assert(theta >= 0 && theta < MAX_ALLOWABLE_HOLE_PAIR_ANGLE);
            const Scalar extra_search_dist = r_hole * tan(theta);
            assert(extra_search_dist > -SMALL_SCALAR);
            return (extra_search_dist + SMALL_SCALAR > -dist) ? 0 : -1;
        }
    } else {
        assert(dist > 0);
        if (ie_h == Ne_hole - 1) {
            return 1;
        } else {
            const Vec3 &x_hole_next = x_hole[ie_h + 2];
            const Vec3 t_next = (x_hole_next - x_hole_B).normalized();
            const Scalar theta = acos(t.dot(t_next));
            assert(theta >= 0 && theta < MAX_ALLOWABLE_HOLE_PAIR_ANGLE);
            const Scalar extra_search_dist = r_hole * tan(theta);
            assert(extra_search_dist > -SMALL_SCALAR);
            return (extra_search_dist + SMALL_SCALAR > dist) ? 0 : 1;
        }
    }
#undef MAX_ALLOWABLE_HOLE_PAIR_ANGLE
}

template <bool initial_search>
inline void update_hole_contact_indices(const Index N, const Index N_hole, const vector<Vec3> &x_hole,
                                        const vector<Scalar> &r_e_hole, vector<Index> &hole_index,
                                        const vector<Vec3> &X, const vector<Vec3> &d_trans) {
    const Index Ne_hole = N_hole - 1;

    if constexpr (initial_search) { /*If its the initial search (all indices are zero from before), the search is more
                                  involved*/
        for (Index i = 0; i < N; i++) {
            int hi = hole_index[i];
            assert(hi == 0); /*Initial value of hi*/

            int search_dir = node_within_hole_segment(hi, N_hole, x_hole, r_e_hole[hi], X[i] + d_trans[i]);

            if (search_dir == 0) {
                continue;
            } else if (search_dir == -1) {
                /*Search backwards*/
                while (search_dir != 0) {
                    assert(search_dir == -1);
                    hi--;
                    assert(hi >= 0);
                    if (hi < 0) {
                        throw runtime_error("Hole index (" + to_string(hi) + ") < num hole elements (" +
                                            to_string(Ne_hole) + ") found for node i=" + to_string(i));
                    }
                    search_dir = node_within_hole_segment(hi, N_hole, x_hole, r_e_hole[hi], X[i] + d_trans[i]);
                }
                hole_index[i] =
                    max(0, hi); /*Not ideal, but more robust than allowing negative int casted to unsigned int */
            } else {
                /*Search forwards*/
                while (search_dir != 0) {
                    assert(search_dir == 1);
                    hi++;
                    assert(hi < (int)Ne_hole);
                    if (hi >= (int)Ne_hole) {
                        throw runtime_error("Hole index (" + to_string(hi) + ") >= num hole elements (" +
                                            to_string(Ne_hole) + ") found for node i=" + to_string(i));
                    }
                    search_dir = node_within_hole_segment(hi, N_hole, x_hole, r_e_hole[hi], X[i] + d_trans[i]);
                }
                hole_index[i] = min(hi, int(Ne_hole - 1));
            }
        }
    } else { /*In later searches its assumed that the hole index is updated on the first try,
           since the node will cross no more than one hole element boundaries per timestep*/
#pragma omp parallel for
        for (Index i = 0; i < N; i++) {
            int hi = hole_index[i];
            assert(hi >= 0 && hi < Ne_hole);
            const int search_dir = node_within_hole_segment(hi, N_hole, x_hole, r_e_hole[hi], X[i] + d_trans[i]);
            hi += search_dir;
            assert(hi >= 0 && hi < Ne_hole);
            hi = max(0, min(hi, int(Ne_hole - 1))); /*Should ideally not happen, but this is done to increase robustness
                                                       and avoid avoid segaults etc. */
            hole_index[i] = hi;
        }
    }
}

inline Scalar dX_e(Index ie, const vector<Vec3> &X) {
    return (X[ie + 1] - X[ie]).norm();
}

inline Scalar dX_element_avg(Index i, const vector<Vec3> &X, const Index N) {
    assert(i < N);
    const Index Ne = N - 1;
    if (i == 0)
        return dX_e(0, X) / 2;
    else if (i == N - 1)
        return dX_e(Ne - 1, X) / 2;
    else
        return (dX_e(i - 1, X) + dX_e(i, X)) / 2;
}

void initialize_hole_contact(const Config &config, const Geometry &geometry, const Borehole &borehole,
                             BeamSystem &beam) {
    if (!config.contact_enabled)
        return;

    constexpr bool initial_search = true;
    try {
        printf("Starting initial hole contact detection. N hole nodes = %i, N string nodes = %i\n",
               borehole.get_N_hole_nodes(), geometry.get_N());

        update_hole_contact_indices<initial_search>(geometry.get_N(), borehole.get_N_hole_nodes(), borehole.get_x(),
                                                    borehole.get_r_hole_element(), beam.hole_index, geometry.get_X(),
                                                    beam.d_trans);
    } catch (const exception &e) {
        throw runtime_error("Hole contact initialization failed\n" + string(e.what()));
    }
}