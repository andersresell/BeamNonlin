// #include "../include/Solver.hpp"

#include "../include/Config.hpp"

void velocity_update_partial(Scalar dt, Index N, const Scalar *__restrict__ M, const Vec3 *__restrict__ J_u,
                             const Vec3 *__restrict__ R_int_trans, const Vec3 *__restrict__ R_int_rot,
                             const Vec3 *__restrict__ R_ext_trans, const Vec3 *__restrict__ R_ext_rot,
                             Vec3 *__restrict__ v_trans, Vec3 *__restrict__ v_rot)
{
#pragma omp parallel
    {
#pragma omp for nowait
        for (Index i = 0; i < N; i++)
        {
            v_trans[i] += 0.5 * dt * (R_ext_trans[i] - R_int_trans[i]) / M[i];

            DEBUG_ONLY(
                if (i == 1) cout << "v_trans:\n"
                                 << v_trans[i] << endl;);
        }
#pragma omp for
        for (Index i = 0; i < N; i++)
        {
            cout << "remove!\n";
            if (i == 0)
                continue;
            /*omega_u_dot = J_u^-1 * (R_rot_u - S(omega_u)*J_u*omega_u)*/
            const Vec3 R_rot_u = R_ext_rot[i] - R_int_rot[i];

            DEBUG_ONLY(if (i == 1) cout
                       << "R_ext_rot\n"
                       << R_ext_rot[i] << endl
                       << "R_int_rot\n"
                       << R_int_rot[i] << endl
                       << "R_rot\n"
                       << R_rot_u << endl);
            Vec3 &omega_u = v_rot[i];

            cout << "R_ext_rot \n"
                 << R_ext_rot[i] << endl;
            // omega_u = {1.424, -2.13, 3.11331};

            Mat3 JJu = J_u[i].asDiagonal();
            Vec3 Ju = J_u[i];
            // Vec3 rot_term_orig = omega_u.cross(Vec3{J_u[i].array() * omega_u.array()});

            // Vec3 rot_term_new = skew_symmetric(omega_u) * JJu * omega_u;

            // cout << "rot_term_orig " << rot_term_orig << endl;
            // cout << "rot_term_new " << rot_term_new << endl;

            cout << "omega_u\n"
                 << omega_u << endl;

            Vec3 omega_u_dot_new;
            omega_u_dot_new.x() = (R_rot_u.x() - (Ju.z() - Ju.y()) * omega_u.y() * omega_u.z()) / Ju.x();
            omega_u_dot_new.y() = (R_rot_u.y() - (Ju.x() - Ju.z()) * omega_u.x() * omega_u.z()) / Ju.y();
            omega_u_dot_new.z() = (R_rot_u.z() - (Ju.y() - Ju.x()) * omega_u.x() * omega_u.y()) / Ju.z();
            // omega_u_dot_new.x() = (R_rot_u.x()) / Ju.x();
            // omega_u_dot_new.y() = (R_rot_u.y()) / Ju.y();
            // omega_u_dot_new.z() = (R_rot_u.z()) / Ju.z();

            // const Vec3 omega_u_dot = (R_rot_u - rot_term_orig).array() / J_u[i].array();
            Vec3 omega_u_dot = omega_u_dot_new;

            cout << "omega_u_dot\n"
                 << omega_u_dot << endl;

            // if (isnan(diff.norm()))
            // {
            //     cout << "nan" << endl;
            //     exit(1);
            // }
            // const Vec3 omega_u_dot = (R_rot_u - omega_u.cross(Vec3{J_u[i].array() * omega_u.array()})).array() / J_u[i].array();
            //  Scalar tol = 0.00001;
            //  if (i != 0 && abs(omega_u[0]) > tol)
            //  {
            //      cout << "omega_u \n"
            //           << omega_u << endl;
            //      assert(false);
            //  }
            Scalar mag_old = omega_u.norm();

            omega_u += 0.5 * dt * omega_u_dot;

            Scalar mag = omega_u.norm();
            Scalar err = mag - mag_old;
            cout << "err " << err << endl;

            cout << "dt " << dt << endl;
            cout << "omega_u updated\n"
                 << omega_u << endl;
            DEBUG_ONLY(
                if (i == 1) cout << "omega_u:\n"
                                 << v_rot[i] << endl;);
            // if (n_glob == 10000 && i > 10)
            //     omega_u[0] = 0.1;
        }
    }
}