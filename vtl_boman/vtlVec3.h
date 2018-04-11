#ifndef VTLVEC3_H
#define VTLVEC3_H

#include "vtl.h"
#include <iostream>

namespace vtl
{

// 
/**
 * @brief A simple vector with 3 components of type T
 */

template <typename T>
class Vec3
{
    T v[3];

  public:
    explicit Vec3(T a1 = 0, T a2 = 0, T a3 = 0)
    {
        v[0] = a1;
        v[1] = a2;
        v[2] = a3;
    }
    Vec3(T a[3])
    {
        v[0] = a[0];
        v[1] = a[1];
        v[2] = a[2];
    }
    friend std::ostream &operator<<(std::ostream &out, Vec3<T> const &obj)
    {
        out << "(" << obj.v[0] << ", " << obj.v[1] << ", " << obj.v[2] << ")";
        return out;
    }
    Vec3<T> operator-(Vec3<T> const &a) const
    {
        return Vec3<T>(v[0] - a.v[0], v[1] - a.v[1], v[2] - a.v[2]);
    }
    Vec3<T> operator+(Vec3<T> const &a) const
    {
        return Vec3<T>(v[0] + a.v[0], v[1] + a.v[1], v[2] + a.v[2]);
    }
    Vec3<T> operator-(T const &a) const
    {
        return Vec3<T>(v[0] - a, v[1] - a, v[2] - a);
    }
    Vec3<T> operator+(T const &a) const
    {
        return Vec3<T>(v[0] + a, v[1] + a, v[2] + a);
    }
    Vec3<T> operator/(T const &d) const
    {
        if (d == 0)
            throw std::runtime_error("division by 0!\n");
        return Vec3<T>(v[0] / d, v[1] / d, v[2] / d);
    }
    Vec3<T> operator/(Vec3<T> const &d) const
    {
        if (d.v[0] == 0 || d.v[1] == 0 || d.v[2] == 0)
            throw std::runtime_error("division by 0!\n");
        return Vec3<T>(v[0] / d.v[0], v[1] / d.v[1], v[2] / d.v[2]);
    }
    Vec3<T> operator*(T const &d) const
    {
        return Vec3<T>(v[0] * d, v[1] * d, v[2] * d);
    }
    Vec3<T> operator*(Vec3<T> const &d) const
    {
        return Vec3<T>(v[0] * d[0], v[1] * d[1], v[2] * d[2]);
    }
    friend Vec3<T> operator*(T const &d, const Vec3<T> &v)
    {
        return v * d;
    }
    Vec3<T> operator-() const
    {
        return Vec3<T>(-v[0], -v[1], -v[2]);
    }
    template <class U>
    operator Vec3<U>(void) const // cast
    {
        return Vec3<U>(v[0], v[1], v[2]);
    }
    T &operator[](int i)
    {
        return v[i];
    }
    T operator[](int i) const
    {
        return v[i];
    }
};

typedef Vec3<double> Vec3d;
typedef Vec3<int> Vec3i;
}

#endif //VTLVEC3_H