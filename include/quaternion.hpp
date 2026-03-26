#ifndef QUATERNION_HPP
#define QUATERNION_HPP

#include "vector3.hpp"
#include <ostream>

struct Quaternion {
    float x, y, z, w;

    Quaternion();

    Quaternion(float x, float y, float z, float w);

    Quaternion operator+(const Quaternion &other) const;

    Quaternion operator-(const Quaternion &other) const;

    Quaternion operator*(const Quaternion &other) const;

    Quaternion operator*(float scalar) const;

    bool operator==(const Quaternion &other) const;

    bool operator!=(const Quaternion &other) const;

    Vector3 operator*(const Vector3 &v) const;

    [[nodiscard]] float length() const;

    [[nodiscard]] float dot(const Quaternion &other) const;

    [[nodiscard]] Quaternion normalized() const;

    [[nodiscard]] Quaternion conjugate() const;

    [[nodiscard]] Quaternion inverse() const;

    [[nodiscard]] Quaternion interpolate(const Quaternion &other, float factor) const;

    [[nodiscard]] bool approximately(const Quaternion &other, float epsilon = 1e-5f) const;

    [[nodiscard]] Vector3 euler() const;

    static Quaternion euler(const Vector3 &degrees);

    static Quaternion between(const Vector3 &from, const Vector3 &to);

    static Quaternion look(const Vector3 &forward, const Vector3 &up = Vector3::UP);

    static Quaternion around(float angle, const Vector3 &axis);

    static const Quaternion IDENTITY;

    friend std::ostream &operator<<(std::ostream &os, const Quaternion &q);
};

#endif