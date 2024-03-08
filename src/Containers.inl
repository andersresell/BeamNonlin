#pragma once
#include "../include/Containers.hpp"

inline Mat3 Quaternion::to_matrix() const
{
    assert(is_close(this->norm(), 1.0));
    Scalar q1 = q.x();
    Scalar q2 = q.y();
    Scalar q3 = q.z();
    Mat3 triad = 2 * Mat3{{q0 * q0 + q1 * q1 - 0.5, q1 * q2 - q3 * q0, q1 * q3 + q2 * q0},
                          {q2 * q1 + q3 * q0, q0 * q0 + q2 * q2 - 0.5, q2 * q3 - q1 * q0},
                          {q3 * q1 - q2 * q0, q3 * q2 + q1 * q0, q0 * q0 + q3 * q3 - 0.5}};

    assert(is_orthogonal(triad));
    return triad;
}

inline void Quaternion::from_matrix(const Mat3 &R)
{

    assert(is_orthogonal(R));
    // See chapter 16.10 crisfield
    Scalar R11 = R(0, 0);
    Scalar R22 = R(1, 1);
    Scalar R33 = R(2, 2);
    Scalar R12 = R(0, 1);
    Scalar R21 = R(1, 0);
    Scalar R23 = R(1, 2);
    Scalar R32 = R(2, 1);
    Scalar R13 = R(0, 2);
    Scalar R31 = R(2, 0);
    Scalar trR = R11 + R22 + R33;
    Scalar &q1 = this->q.x();
    Scalar &q2 = this->q.y();
    Scalar &q3 = this->q.z();

    Scalar a = max(max(trR, R11), max(R22, R33));
    if (a == trR)
    {
        q0 = 0.5 * sqrt(1 + a);
        q1 = (R32 - R23) / (4 * q0);
        q2 = (R13 - R31) / (4 * q0);
        q3 = (R21 - R12) / (4 * q0);
    }
    else if (a == R11)
    {
        q1 = sqrt(0.5 * a + 0.25 * (1 - trR));
        q0 = 0.25 * (R32 - R23) / q1;
        q2 = 0.25 * (R21 + R12) / q1;
        q3 = 0.25 * (R31 + R13) / q1;
    }
    else if (a == R22)
    {
        q2 = sqrt(0.5 * a + 0.25 * (1 - trR));
        q0 = 0.25 * (R13 - R31) / q2;
        q1 = 0.25 * (R12 + R21) / q2;
        q3 = 0.25 * (R32 + R23) / q2;
    }
    else
    {
        assert(a == R33);
        q3 = sqrt(0.5 * a + 0.25 * (1 - trR));
        q0 = 0.25 * (R21 - R12) / q3;
        q1 = 0.25 * (R13 + R31) / q3;
        q2 = 0.25 * (R23 + R32) / q3;
    }
    assert(is_close(norm(), 1.0));
    this->normalize(); // needed ?
}

inline Scalar Quaternion::norm() const
{
    return sqrt(norm_sqr());
}
inline Scalar Quaternion::norm_sqr() const
{
    return q0 * q0 + q.x() * q.x() + q.y() * q.y() + q.z() * q.z();
}

inline void Quaternion::normalize()
{
    Scalar inv_norm = 1.0 / this->norm();
    q0 *= inv_norm; // probably not faster than just divding by norm after compiler optimizations
    q.x() *= inv_norm;
    q.y() *= inv_norm;
    q.z() *= inv_norm;
}

inline Quaternion::Quaternion(const Vec3 &Theta)
{
    Scalar theta = Theta.norm();
    q0 = cos(theta / 2);
    q = sin(theta / 2) * Theta.normalized();
    assert(is_close(this->norm_sqr(), 1));
}

inline Quaternion Quaternion::product(const Quaternion &a) const
{
    assert(is_close(a.norm_sqr(), 1));
    /*16.72 crisfield*/
    const Quaternion &b{*this};
    Quaternion q_ab;
    q_ab.q0 = a.q0 * b.q0 - a.q.dot(b.q);
    q_ab.q = a.q0 * b.q + b.q0 * a.q - a.q.cross(b.q);
    assert(is_close(q_ab.norm_sqr(), 1));
    return q_ab;
}

inline void Quaternion::compound_rotate(const Vec3 &Theta)
{
    Quaternion delta_q{Theta};
    *this = delta_q.product(*this);
}

inline Vec3 Quaternion::rotate_vector(const Vec3 &v0) const
{ /*Simply multiplying Rodrigues formula for quaternions with the vector*/
    const Scalar q1 = q.x();
    const Scalar q2 = q.y();
    const Scalar q3 = q.z();
    Vec3 vn;
    vn.x() = 2 * ((q0 * q0 + q1 * q1 - 0.5) * v0.x() + (q1 * q2 - q3 * q0) * v0.y() + (q1 * q3 + q2 * q0) * v0.z());
    vn.y() = 2 * ((q2 * q1 + q3 * q0) * v0.x() + (q0 * q0 + q2 * q2 - 0.5) * v0.y() + (q2 * q3 - q1 * q0) * v0.z());
    vn.z() = 2 * ((q3 * q1 - q2 * q0) * v0.x() + (q3 * q2 + q1 * q0) * v0.y() + (q0 * q0 + q3 * q3 - 0.5) * v0.z());
    assert(is_close(vn.norm(), v0.norm()));
    return vn;
}

inline Vec3 Quaternion::rotate_vector_reversed(const Vec3 &v0) const
{ /*Simply multiplying the transposed of Rodrigues formula for quaternions with the vector*/
    const Scalar q1 = q.x();
    const Scalar q2 = q.y();
    const Scalar q3 = q.z();
    Vec3 vn;
    vn.x() = 2 * ((q0 * q0 + q1 * q1 - 0.5) * v0.x() + (q2 * q1 + q3 * q0) * v0.y() + (q3 * q1 - q2 * q0) * v0.z());
    vn.y() = 2 * ((q1 * q2 - q3 * q0) * v0.x() + (q0 * q0 + q2 * q2 - 0.5) * v0.y() + (q3 * q2 + q1 * q0) * v0.z());
    vn.z() = 2 * ((q1 * q3 + q2 * q0) * v0.x() + (q2 * q3 - q1 * q0) * v0.y() + (q0 * q0 + q3 * q3 - 0.5) * v0.z());
    assert(is_close(this->norm_sqr(), 1));
    assert(is_close(vn.norm(), v0.norm()));
    return vn;
}
// inline Vec3Quat::Vec3Quat()
// {
//     static_assert(sizeof(Vec3Quat) == sizeof(Scalar) * 7);
//     trans = {0, 0, 0};
//     rot.from_matrix(Mat3::Identity());
// }

// inline void Vec3Quat::print_array(vector<Vec3Quat> arr, string label, bool print_trans, bool print_rot)
// {
//     if (label != "")
//         cout << label << ":" << endl;
//     if (print_trans)
//     {
//         cout << "trans:\n";
//         for (const auto &e : arr)
//         {
//             cout << e.trans.transpose() << ",\n";
//         }
//     }
//     if (print_rot)
//     {
//         cout << "\nrot:\n";
//         for (const auto &e : arr)
//         {
//             cout << "{" << e.rot.to_matrix() << "},\n";
//         }
//     }
//     cout << endl;
// }

// inline ostream &operator<<(ostream &os, const Vec3Quat &rhs)
// {
//     return os << "{" << rhs.trans.transpose() << ",\n"
//               << rhs.rot.to_matrix() << "}";
// }