#ifndef POINT_ALIGN_H
#define POINT_ALIGN_H

#include <math.h>

#include "hashmap.h"
#define PI 3.14159265358979323846

typedef enum {
    RETURN_OK,
    RETURN_FAIL,
    RETURN_INVERSE_FAIL,
    RETURN_ITERATION_FAIL,
} ReturnCode;

typedef struct Vec3d {
    double x, y, z;
} Vec3d;
#define VEC3D_ZERO {0.0, 0.0, 0.0}

typedef struct Vec4d {
    double x, y, z, w;
} Vec4d;
#define VEC4D_ZERO {0.0, 0.0, 0.0, 0.0}

typedef struct Vec4d Quatd;
#define QUATD_IDENTITY {0.0, 0.0, 0.0, 1.0}

/* 4x4 matrix with column-major memory order */
typedef struct Mat4d {
    double xi, xj, xk, xw;
    double yi, yj, yk, yw;
    double zi, zj, zk, zw;
    double ti, tj, tk, tw;
} Mat4d;
// clang-format off
#define MAT4D_IDENTITY { \
    1.0, 0.0, 0.0, 0.0, \
    0.0, 1.0, 0.0, 0.0, \
    0.0, 0.0, 1.0, 0.0, \
    0.0, 0.0, 0.0, 1.0, \
}
// clang-format on

typedef struct FerrariCoeffs {
    double values[5];
} FerrariCoeffs;

typedef struct Transform3d {
    Vec3d translation;
    Quatd rotation;
} Transform3d;

DEFINE_HASHMAP(Vec3d)

Vec3d compute_centroid(const Vec3d* points, uint32_t count);
Vec3d transform_point(const Vec3d pt, const Transform3d t);
double compute_rmsd(
    const Vec3d* source, const Vec3d* target, uint32_t n_points, Transform3d t
);
ReturnCode Mat4d_inverse(const Mat4d m, Mat4d* dest);
Transform3d align_points(
    const Vec3d* source, const Vec3d* target, uint32_t n_points
);

/* INLINE FUNCTIONS */

static inline Vec3d Vec3d_init(double x, double y, double z) {
    Vec3d v = {x, y, z};
    return v;
}
static inline Vec3d Vec3d_add(Vec3d va, Vec3d vb) {
    Vec3d result = {va.x + vb.x, va.y + vb.y, va.z + vb.z};
    return result;
}

static inline Vec3d Vec3d_sub(Vec3d va, Vec3d vb) {
    Vec3d result = {va.x - vb.x, va.y - vb.y, va.z - vb.z};
    return result;
}

static inline Vec3d Vec3d_scale(Vec3d va, double s) {
    Vec3d result = {va.x * s, va.y * s, va.z * s};
    return result;
}

static inline double Vec3d_dot(Vec3d va, Vec3d vb) {
    return va.x * vb.x + va.y * vb.y + va.z * vb.z;
}

static inline double Vec3d_length(Vec3d v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

static inline Vec3d Vec3d_cross(Vec3d a, Vec3d b) {
    Vec3d c = VEC3D_ZERO;
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    return c;
}

static inline double Quatd_dot(Quatd qa, Quatd qb) {
    return qa.x * qb.x + qa.y * qb.y + qa.z * qb.z + qa.w * qb.w;
}

static inline Quatd Quatd_mul(Quatd qa, Quatd qb) {
    Quatd result = QUATD_IDENTITY;

    result.x = qa.x * qb.w + qa.w * qb.x + qa.y * qb.z - qa.z * qb.y;
    result.y = qa.y * qb.w + qa.w * qb.y + qa.z * qb.x - qa.x * qb.z;
    result.z = qa.z * qb.w + qa.w * qb.z + qa.x * qb.y - qa.y * qb.x;
    result.w = qa.w * qb.w - qa.x * qb.x - qa.y * qb.y - qa.z * qb.z;

    return result;
}

static inline double Quatd_length(Quatd q) {
    return sqrt(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
}

static inline Quatd Quatd_normalize(Quatd q) {
    double mag = Quatd_length(q);
    if (mag == 0.0) {
        mag = 1.0;
    }
    double imag = 1.0 / mag;
    q.x = q.x * imag;
    q.y = q.y * imag;
    q.z = q.z * imag;
    q.w = q.w * imag;
    return q;
}

static inline Vec3d Quatd_rotate_vec3d(Quatd q, Vec3d v) {
    // Quaternion rotation: v' = v + 2 * cross(q.xyz, cross(q.xyz, v) + q.w * v)
    Vec3d qv = Vec3d_init(q.x, q.y, q.z);
    Vec3d t = Vec3d_cross(qv, v);
    t = Vec3d_add(t, Vec3d_scale(v, q.w));
    t = Vec3d_cross(qv, t);
    t = Vec3d_scale(t, 2.0);
    return Vec3d_add(v, t);
}

/** Pre-multiply a `Vec4d` by a `Mat4d`: result = M * v
 *
 * Note that this makes an independent copy of `m` and `v`
 */
static inline Vec4d Mat4d_mul_vec4d(Mat4d m, Vec4d v) {
    Vec4d result;
    result.x = m.xi * v.x + m.yi * v.y + m.zi * v.z + m.ti * v.w;
    result.y = m.xj * v.x + m.yj * v.y + m.zj * v.z + m.tj * v.w;
    result.z = m.xk * v.x + m.yk * v.y + m.zk * v.z + m.tk * v.w;
    result.w = m.xw * v.x + m.yw * v.y + m.zw * v.z + m.tw * v.w;
    return result;
}

static inline double det2(double a, double b, double c, double d) {
    return a * d - b * c;
}

static inline double det3(
    double a,
    double b,
    double c,
    double d,
    double e,
    double f,
    double g,
    double h,
    double i
) {
    return a * det2(e, f, h, i) - b * det2(d, f, g, i) + c * det2(d, e, g, h);
}
static inline Mat4d Mat4d_sub_scalar_diag(Mat4d m, double lambda) {
    // Subtract lambda from diagonal elements only
    m.xi -= lambda;
    m.yj -= lambda;
    m.zk -= lambda;
    m.tw -= lambda;

    return m;
}

#endif /* POINT_ALIGN_H */
