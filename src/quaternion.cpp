#include "quaternion.hpp"
#include <algorithm>

const Quaternion Quaternion::IDENTITY = {0, 0, 0, 1};

Quaternion::Quaternion() : x(0), y(0), z(0), w(1) {}

Quaternion::Quaternion(const float x, const float y, const float z, const float w) : x(x), y(y), z(z), w(w) {}

Quaternion Quaternion::operator+(const Quaternion &other) const {
    return {x + other.x, y + other.y, z + other.z, w + other.w};
}

Quaternion Quaternion::operator-(const Quaternion &other) const {
    return {x - other.x, y - other.y, z - other.z, w - other.w};
}

Quaternion Quaternion::operator*(const Quaternion &other) const {
    return {
        w * other.x + x * other.w + y * other.z - z * other.y,
        w * other.y - x * other.z + y * other.w + z * other.x,
        w * other.z + x * other.y - y * other.x + z * other.w,
        w * other.w - x * other.x - y * other.y - z * other.z
    };
}

Quaternion Quaternion::operator*(const float scalar) const {
    return {x * scalar, y * scalar, z * scalar, w * scalar};
}

bool Quaternion::operator==(const Quaternion &other) const {
    return x == other.x && y == other.y && z == other.z && w == other.w;
}

bool Quaternion::operator!=(const Quaternion &other) const {
    return !(*this == other);
}

Vector3 Quaternion::operator*(const Vector3 &v) const {
    const Vector3 q(x, y, z);
    const Vector3 crossed = q.cross(v);
    return v + crossed * (2.0f * w) + q.cross(crossed) * 2.0f;
}

float Quaternion::length() const {
    return sqrtf(x * x + y * y + z * z + w * w);
}

float Quaternion::dot(const Quaternion &other) const {
    return x * other.x + y * other.y + z * other.z + w * other.w;
}

Quaternion Quaternion::normalized() const {
    const float len = length();
    if (len == 0) return IDENTITY;
    return {x / len, y / len, z / len, w / len};
}

Quaternion Quaternion::conjugate() const {
    return {-x, -y, -z, w};
}

Quaternion Quaternion::inverse() const {
    const float magnitude = dot(*this);
    if (magnitude == 0) return IDENTITY;
    return {-x / magnitude, -y / magnitude, -z / magnitude, w / magnitude};
}

Quaternion Quaternion::interpolate(const Quaternion &other, float factor) const {
    factor = std::clamp(factor, 0.0f, 1.0f);
    float cosine = dot(other);
    Quaternion target = other;

    if (cosine < 0.0f) {
        cosine = -cosine;
        target = {-other.x, -other.y, -other.z, -other.w};
    }

    if (cosine > 0.9995f) {
        return (*this + (target - *this) * factor).normalized();
    }

    const float base  = acosf(cosine);
    const float scale = sinf(base);
    return (*this * sinf(base * (1.0f - factor)) + target * sinf(base * factor)) * (1.0f / scale);
}

bool Quaternion::approximately(const Quaternion &other, const float epsilon) const {
    return fabsf(dot(other)) >= 1.0f - epsilon;
}

Vector3 Quaternion::euler() const {
    Vector3 angles;

    const float sinX = 2.0f * (w * x + y * z);
    const float cosX = 1.0f - 2.0f * (x * x + y * y);
    angles.x = atan2f(sinX, cosX) * (180.0f / 3.14159265f);

    const float sinY = 2.0f * (w * y - z * x);
    angles.y = fabsf(sinY) >= 1.0f
        ? copysignf(90.0f, sinY)
        : asinf(sinY) * (180.0f / 3.14159265f);

    const float sinZ = 2.0f * (w * z + x * y);
    const float cosZ = 1.0f - 2.0f * (y * y + z * z);
    angles.z = atan2f(sinZ, cosZ) * (180.0f / 3.14159265f);

    return angles;
}

Quaternion Quaternion::euler(const Vector3 &degrees) {
    const float rad = 3.14159265f / 180.0f;
    const float hx  = degrees.x * rad * 0.5f;
    const float hy  = degrees.y * rad * 0.5f;
    const float hz  = degrees.z * rad * 0.5f;

    const float sx = sinf(hx), cx = cosf(hx);
    const float sy = sinf(hy), cy = cosf(hy);
    const float sz = sinf(hz), cz = cosf(hz);

    return {
        cx * sy * sz + sx * cy * cz,
        cx * sy * cz - sx * cy * sz,
        cx * cy * sz + sx * sy * cz,
        cx * cy * cz - sx * sy * sz
    };
}

Quaternion Quaternion::between(const Vector3 &from, const Vector3 &to) {
    const Vector3 f = from.normalized();
    const Vector3 t = to.normalized();
    const float cosine = f.dot(t);

    if (cosine >= 1.0f - 1e-6f) return IDENTITY;

    if (cosine <= -1.0f + 1e-6f) {
        Vector3 axis = Vector3::RIGHT.cross(f);
        if (axis.length() < 1e-6f) axis = Vector3::UP.cross(f);
        return around(180.0f, axis.normalized());
    }

    const Vector3 axis = f.cross(t).normalized();
    const float angle  = acosf(cosine) * (180.0f / 3.14159265f);
    return around(angle, axis);
}

Quaternion Quaternion::look(const Vector3 &forward, const Vector3 &up) {
    const Vector3 f = forward.normalized();
    const Vector3 r = up.cross(f).normalized();
    const Vector3 u = f.cross(r);

    const float trace = r.x + u.y + f.z;

    if (trace > 0.0f) {
        const float root = sqrtf(trace + 1.0f) * 2.0f;
        return {
            (u.z - f.y) / root,
            (f.x - r.z) / root,
            (r.y - u.x) / root,
            root * 0.25f
        };
    }

    if (r.x > u.y && r.x > f.z) {
        const float root = sqrtf(1.0f + r.x - u.y - f.z) * 2.0f;
        return { root * 0.25f, (r.y + u.x) / root, (f.x + r.z) / root, (u.z - f.y) / root };
    }

    if (u.y > f.z) {
        const float root = sqrtf(1.0f + u.y - r.x - f.z) * 2.0f;
        return { (r.y + u.x) / root, root * 0.25f, (u.z + f.y) / root, (f.x - r.z) / root };
    }

    const float root = sqrtf(1.0f + f.z - r.x - u.y) * 2.0f;
    return { (f.x + r.z) / root, (u.z + f.y) / root, root * 0.25f, (r.y - u.x) / root };
}

Quaternion Quaternion::around(const float angle, const Vector3 &axis) {
    const float rad  = angle * (3.14159265f / 180.0f) * 0.5f;
    const float sinA = sinf(rad);
    const Vector3 a  = axis.normalized();
    return {a.x * sinA, a.y * sinA, a.z * sinA, cosf(rad)};
}

std::ostream &operator<<(std::ostream &os, const Quaternion &q) {
    os << "(" << q.x << ", " << q.y << ", " << q.z << ", " << q.w << ")";
    return os;
}