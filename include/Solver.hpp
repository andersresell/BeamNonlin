#pragma once
#include "Config.hpp"
#include "Containers.hpp"
#include "Geometry.hpp"
#include "SolverUtils.hpp"

void solve(Config &config, const Geometry &geometry);

void calc_dt(Config &config, const Geometry &geometry);

inline void check_energy_balance(const Config &config, const BeamSystem &beam_sys);

static Index n_glob;

inline void assemble(const Config &config, const Geometry &geometry, BeamSystem &beam_sys);

inline void calc_element_inner_forces(Index ie, const Vec3 *__restrict__ X, const Vec3 *__restrict__ d_trans,
                                      const Quaternion *__restrict__ d_rot,
                                      Vec3Vec3 *__restrict__ R_int, Scalar ri_e, Scalar ro_e, Scalar E, Scalar G);

inline Vec7 calc_nodal_forces_local(Scalar ri, Scalar ro, Scalar l0, Scalar E, Scalar G, Scalar ul,
                                    Scalar theta_1l, Scalar theta_2l, Scalar theta_3l, Scalar theta_4l,
                                    Scalar theta_5l, Scalar theta_6l);

// inline void step_central_differences(Scalar dt, Index N, Vec3Quat *__restrict__ u, Vec3Vec3 *__restrict__ v,
//                                      const Scalar *__restrict__ M, const Vec3 *__restrict__ J_u,
//                                      const Vec3Vec3 *__restrict__ R_int, const Vec3Vec3 *__restrict__ R_ext,
//                                      bool check_energy_balance, Scalar &W_int, Scalar &W_ext, Scalar &KE);

inline void calc_static_loads(const Config &config, const Geometry &geometry, vector<Vec3Vec3> &R_static);

inline void set_simple_bc(BC_Case bc_case, const Geometry &geometry, BeamSystem &beam_system);

inline void velocity_update_partial(Scalar dt, Index N, const Scalar *__restrict__ M,
                                    const Vec3 *__restrict__ J_u, const Vec3Vec3 *__restrict__ R_int,
                                    const Vec3Vec3 *__restrict__ R_ext,
                                    Vec3Vec3 *__restrict__ v)
{
#pragma omp parallel for
    for (Index i = 0; i < N; i++)
    {
        v[i].trans += 0.5 * dt * (R_ext[i].trans - R_int[i].trans) / M[i];

        /*omega_u_dot = J_u^-1 * (R_rot_u - S(omega_u)*J_u*omega_u)*/
        const Vec3 R_rot_u = R_ext[i].rot - R_int[i].rot;
        const Vec3 &omega_u = v[i].rot;
        const Vec3 omega_u_dot = (R_rot_u - omega_u.cross(Vec3{J_u[i].array() * omega_u.array()})).array() / J_u[i].array();
        v[i].rot += 0.5 * dt * omega_u_dot;
    }
}

inline void displacement_update(Scalar dt, Index N, Vec3Vec3 *__restrict__ v, Vec3 *__restrict__ d_trans,
                                Quaternion *__restrict__ d_rot)
{
#pragma omp parallel for
    for (Index i = 0; i < N; i++)
    {
        d_trans[i] += dt * v[i].trans;
        Mat3 U = d_rot[i].to_matrix();
        U = U * (Mat3::Identity() - 0.5 * dt * skew_symmetric(v[i].rot)).inverse() *
            (Mat3::Identity() + 0.5 * dt * skew_symmetric(v[i].rot));
        assert(U.allFinite());
        d_rot[i].from_matrix(U);
    }
}

inline void calc_delta_d(Scalar dt, Index N, Vec3Vec3 *__restrict__ delta_d, const Vec3Vec3 *__restrict__ v)
{
#pragma omp parallel for
    for (Index i = 0; i < N; i++)
    {
        delta_d[i].trans = dt * v[i].trans;
        delta_d[i].rot = dt * v[i].rot;
    }
}

inline void work_update_partial(Index N, const Vec3Vec3 *__restrict__ delta_d,
                                const Vec3Vec3 *__restrict__ R_int, const Vec3Vec3 *__restrict__ R_ext,
                                Scalar &W_ext, Scalar &W_int)
{
#pragma omp parallel for reduction(+ : W_int) reduction(+ : W_ext)
    for (Index i = 0; i < N; i++)
    {
        /*The rotational dofs should have been rotated to the body frame allready*/
        W_int += 0.5 * (delta_d[i].trans.dot(R_int[i].trans) + delta_d[i].rot.dot(R_int[i].rot));
        W_ext += 0.5 * (delta_d[i].trans.dot(R_ext[i].trans) + delta_d[i].rot.dot(R_ext[i].rot));
    }
}

inline void update_kinetic_energy(Index N, const Scalar *__restrict__ M, const Vec3 *__restrict__ J_u,
                                  const Vec3Vec3 *__restrict__ v, Scalar &KE)
{
    KE = 0;
#pragma omp parallel for reduction(+ : KE)
    for (Index i = 0; i < N; i++)
    {
        const Vec3 &vt = v[i].trans;
        const Vec3 &omega_u = v[i].rot;
        KE += 0.5 * M[i] * vt.squaredNorm() + 0.5 * omega_u.dot(Vec3{J_u[i].array() * omega_u.array()});
    }
}

inline void rotate_R_rot_to_body_frame(Index N, const Quaternion *__restrict__ d_rot,
                                       Vec3Vec3 *__restrict__ R_int, Vec3Vec3 *__restrict__ R_ext)
{
    for (Index i = 0; i < N; i++)
    {
        Mat3 U = d_rot[i].to_matrix();
        R_int[i].rot = U.transpose() * R_int[i].rot;
        R_ext[i].rot = U.transpose() * R_ext[i].rot;
    }
}

inline void step_explicit(Config &config, const Geometry &geometry, BeamSystem &beam_sys)
{
    const Scalar dt = config.dt;
    const Index N = geometry.get_N();
    vector<Vec3> &d_trans = beam_sys.d_trans;
    vector<Quaternion> &d_rot = beam_sys.d_rot;
    vector<Vec3Vec3> &v = beam_sys.v;
    vector<Vec3Vec3> &R_int = beam_sys.R_int;
    vector<Vec3Vec3> &R_ext = beam_sys.R_ext;
    const vector<Scalar> &M = beam_sys.M;
    const vector<Vec3> &J_u = beam_sys.J_u;
    const bool check_energy_balance = config.check_energy_balance;
    Scalar &W_int = beam_sys.W_int;
    Scalar &W_ext = beam_sys.W_ext;
    Scalar &KE = beam_sys.KE;
    vector<Vec3Vec3> &delta_d = beam_sys.delta_d; /*Only used if energy balance is checked*/
    if (check_energy_balance)
    {
        assert(delta_d.size() == N);
    }
    else
    {
        assert(delta_d.size() == 0);
    }

    /*velocity at t_{n+1/2}*/
    velocity_update_partial(dt, N, M.data(), J_u.data(), R_int.data(), R_ext.data(), v.data());

    /*Enforcing boundary conditions*/
    set_simple_bc(config.bc_case, geometry, beam_sys);

    /*Checking energy balance in two steps.
    The equation to be computed is W_{n+1} = W_n + delta_u/2*(Rn + R_{n+1})
    This will be done in two steps in order to not need double storage for R
    First W += delta_u/2*R_n, then updating R, then W +=  * delta_u/2*R_{n+1}
    */
    if (check_energy_balance)
    {
        calc_delta_d(dt, N, delta_d.data(), v.data());

        work_update_partial(N, delta_d.data(), R_int.data(), R_ext.data(), W_ext, W_int);
    }

    /*Displacement at t_{n+1}*/
    displacement_update(dt, N, v.data(), d_trans.data(), d_rot.data());

    /*Updating external and internal forces*/
    assemble(config, geometry, beam_sys);

    /*Since the moments are allways used in the body frame, these are rotated
    once and for all instead of every time they are needed.*/
    rotate_R_rot_to_body_frame(N, d_rot.data(), R_int.data(), R_ext.data());

    if (check_energy_balance)
    {
        work_update_partial(N, delta_d.data(), R_int.data(), R_ext.data(), W_ext, W_int);
    }

    /*velocity at t_{n+1} */
    velocity_update_partial(dt, N, M.data(), J_u.data(), R_int.data(), R_ext.data(), v.data());

    /*Enforcing boundary conditions*/
    set_simple_bc(config.bc_case, geometry, beam_sys);

    if (check_energy_balance)
    {
        update_kinetic_energy(N, M.data(), J_u.data(), v.data(), KE);
    }
}

#include "../src/Solver.inl"
