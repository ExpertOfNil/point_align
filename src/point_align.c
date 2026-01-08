#include "point_align.h"

#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_ITERATIONS 10000
#define CONVERGENCE_TOLERANCE 0.999999999

static void Mat4d_print(Mat4d m) {
    printf("Col 0: [%12.6f,%12.6f,%12.6f,%12.6f]\n", m.xi, m.xj, m.xk, m.xw);
    printf("Col 1: [%12.6f,%12.6f,%12.6f,%12.6f]\n", m.yi, m.yj, m.yk, m.yw);
    printf("Col 2: [%12.6f,%12.6f,%12.6f,%12.6f]\n", m.zi, m.zj, m.zk, m.zw);
    printf("Col 3: [%12.6f,%12.6f,%12.6f,%12.6f]\n", m.ti, m.tj, m.tk, m.tw);
    printf("Diag : [%12.6f,%12.6f,%12.6f,%12.6f]\n", m.xi, m.yj, m.zk, m.tw);
}

ReturnCode Mat4d_inverse(const Mat4d m, Mat4d* dest) {
    Vec3d a = {m.xi, m.xj, m.xk};
    double x = m.xw;
    Vec3d b = {m.yi, m.yj, m.yk};
    double y = m.yw;
    Vec3d c = {m.zi, m.zj, m.zk};
    double z = m.zw;
    Vec3d d = {m.ti, m.tj, m.tk};
    double w = m.tw;

    Vec3d s = Vec3d_cross(a, b);
    Vec3d t = Vec3d_cross(c, d);
    Vec3d u = Vec3d_sub(Vec3d_scale(a, y), Vec3d_scale(b, x));
    Vec3d v = Vec3d_sub(Vec3d_scale(c, w), Vec3d_scale(d, z));

    double det = (Vec3d_dot(s, v) + Vec3d_dot(t, u));
    if (fabs(det) < 1e-12) {
        return RETURN_INVERSE_FAIL;
    }
    double inv_det = 1.0 / det;
    s = Vec3d_scale(s, inv_det);
    t = Vec3d_scale(t, inv_det);
    u = Vec3d_scale(u, inv_det);
    v = Vec3d_scale(v, inv_det);

    Vec3d r0 = Vec3d_add(Vec3d_cross(b, v), Vec3d_scale(t, x));
    Vec3d r1 = Vec3d_sub(Vec3d_cross(v, a), Vec3d_scale(t, y));
    Vec3d r2 = Vec3d_add(Vec3d_cross(d, u), Vec3d_scale(s, w));
    Vec3d r3 = Vec3d_sub(Vec3d_cross(u, c), Vec3d_scale(s, z));

    *dest = (Mat4d){
        // clang-format off
        r0.x, r1.x, r2.x, r3.x,
        r0.y, r1.y, r2.y, r3.y,
        r0.z, r1.z, r2.z, r3.z,
        -Vec3d_dot(b, t), Vec3d_dot(a, t), -Vec3d_dot(d, s), Vec3d_dot(c, s),
        // clang-format on
    };
    return RETURN_OK;
}

/* Compute the centroid of a point array */
Vec3d compute_centroid(const Vec3d* points, uint32_t count) {
    Vec3d centroid = VEC3D_ZERO;
    for (uint32_t i = 0; i < count; ++i) {
        centroid = Vec3d_add(centroid, points[i]);
    }
    return Vec3d_scale(centroid, 1.0 / (double)count);
}

/* Creates a transformed copy of a `Vec3` */
Vec3d transform_point(const Vec3d pt, const Transform3d t) {
    Vec3d rotated = Quatd_rotate_vec3d(t.rotation, pt);
    return Vec3d_add(rotated, t.translation);
}

/* Compute the root-mean-squared deviation */
double compute_rmsd(
    const Vec3d* source, const Vec3d* target, uint32_t n_points, Transform3d t
) {
    if (n_points == 0) return 0.0;
    double sum_sq_dist = 0.0;
    for (size_t i = 0; i < n_points; i++) {
        Vec3d transformed = transform_point(source[i], t);
        Vec3d diff = Vec3d_sub(transformed, target[i]);
        sum_sq_dist += Vec3d_dot(diff, diff);
    }
    return sqrt(sum_sq_dist / (double)n_points);
}

static ReturnCode power_iteration(const Mat4d N, uint32_t max_iters, Quatd* q) {
    Vec4d q_vec = *q;
    uint32_t iter = 0;
    double q_delta = DBL_MAX;
    Quatd convergence[MAX_ITERATIONS] = {0};
    for (; iter < max_iters; iter++) {
        // q_new = N * q
        Vec4d q_new = Mat4d_mul_vec4d(N, q_vec);

        // Normalize
        q_new = Quatd_normalize(q_new);
        convergence[iter] = q_new;

        // Check convergence
        q_delta = Quatd_dot(q_vec, q_new);
        if (q_delta < 0.0) {
            q_new.x = -q_new.x;
            q_new.y = -q_new.y;
            q_new.z = -q_new.z;
            q_new.w = -q_new.w;
            q_delta = -q_delta;
        }
        q_vec = q_new;
        if (q_delta > CONVERGENCE_TOLERANCE) {
            printf("Convergence reached!\n");
            break;  // Converged early
        }
    }
    printf("Power Iterations: %d\n", iter);
    printf(
        "q_vec (x,y,z,w): [%12.6f, %12.6f, %12.6f, %12.6f]\n",
        q_vec.x,
        q_vec.y,
        q_vec.z,
        q_vec.w
    );

    double eigenvalue = Quatd_dot(Mat4d_mul_vec4d(N, q_vec), q_vec);
    printf("Power Eigenvalue: %12.6f\n", eigenvalue);
    if (iter >= MAX_ITERATIONS) {
        FILE* fp = fopen("power_iterations.txt", "w");
        assert(fp != NULL && "Failed to open iterations log\n");
        fprintf(fp, "x, y, z, w\n");
        for (uint32_t i = 0; i < MAX_ITERATIONS; ++i) {
            Quatd q_data = convergence[i];
            fprintf(
                fp,
                "[%12.6f,%12.6f,%12.6f,%12.6f]\n",
                q_data.x,
                q_data.y,
                q_data.z,
                q_data.w
            );
        }
        return RETURN_ITERATION_FAIL;
    }
    *q = q_vec;
    return RETURN_OK;
}

// Solve general cubic: t³ + a*t² + b*t + c = 0
// Returns number of real roots
static int solve_cubic(double a, double b, double c, double roots[3]) {
    // Depress the cubic: substitute t = y - a/3
    // This eliminates the t² term
    double p = b - a * a / 3.0;
    double q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;

    /* Solve depressed cubic: t³ + p*t + q = 0 */
    int num_roots = 0;

    // Compute discriminant
    double discriminant = q * q / 4.0 + p * p * p / 27.0;
    printf("  Cubic discriminant: %.12f\n", discriminant);
    if (discriminant > 1e-12) {
        // One real root (discriminant > 0)
        double sqrt_disc = sqrt(discriminant);
        double u = cbrt(-q / 2.0 + sqrt_disc);
        double v = cbrt(-q / 2.0 - sqrt_disc);
        roots[0] = u + v;
        printf("  Method: One real root (Cardano)\n");
        num_roots = 1;
    } else if (discriminant > -1e-12) {
        // All real, at least two equal (discriminant ≈ 0)
        double u = cbrt(-q / 2.0);
        roots[0] = 2.0 * u;
        roots[1] = -u;
        roots[2] = -u;
        num_roots = 3;
        printf("  Method: All real, at least two equal\n");
    } else {
        // Three distinct real roots (discriminant < 0)
        // Use trigonometric method to avoid complex arithmetic
        double r = 2.0 * sqrt(-p / 3.0);
        double phi = acos(3.0 * q / (p * r));

        roots[0] = r * cos(phi / 3.0);
        roots[1] = r * cos((phi + 2.0 * PI) / 3.0);
        roots[2] = r * cos((phi + 4.0 * PI) / 3.0);
        num_roots = 3;
        printf("  Method: Three distinct real (Trigonometric)\n");
        // After computing r and phi (if you get to this branch):
        printf("  Trig method: r=%.12f, phi=%.12f\n", r, phi);
    }

    /* Transform back: t = y - a/3 */
    double offset = a / 3.0;
    for (int i = 0; i < num_roots; i++) {
        roots[i] -= offset;
    }

    return num_roots;
}

// Solve quartic: λ⁴ + a*λ³ + b*λ² + c*λ + d = 0
// Returns number of real roots (0 to 4)
// Roots stored in roots[0..3], sorted in descending order
static uint32_t solve_quartic(
    double a, double b, double c, double d, double roots[4]
) {
    // Depress the quartic by substituting λ = y - a/4
    double a2 = a * a;
    double a3 = a2 * a;
    double a4 = a2 * a2;

    double p = b - 3.0 * a2 / 8.0;
    double q = c + a3 / 8.0 - a * b / 2.0;
    double r = d - 3.0 * a4 / 256.0 + a2 * b / 16.0 - a * c / 4.0;
    printf(
        "DEBUG: Depressed quartic coefficients: p=%.6f, q=%.6f, r=%.6f\n",
        p,
        q,
        r
    );

    // Special case: depressed quartic is biquadratic (q = 0)
    if (fabs(q) < 1e-12) {
        // Solve y⁴ + py² + r = 0 as quadratic in y²
        double discriminant = p * p - 4.0 * r;
        if (discriminant < 0) return 0;  // No real roots

        double sqrt_disc = sqrt(discriminant);
        double y2_1 = (-p + sqrt_disc) / 2.0;
        double y2_2 = (-p - sqrt_disc) / 2.0;

        int count = 0;
        if (y2_1 >= 0) {
            roots[count++] = sqrt(y2_1) - a / 4.0;
            roots[count++] = -sqrt(y2_1) - a / 4.0;
        }
        if (y2_2 >= 0 && fabs(y2_1 - y2_2) > 1e-12) {
            roots[count++] = sqrt(y2_2) - a / 4.0;
            roots[count++] = -sqrt(y2_2) - a / 4.0;
        }
        return count;
    }

    // Solve resolvent cubic for parameter u
    //
    // The depressed quartic y⁴ + py² + qy + r = 0 can be factored into
    // two quadratics if we find a parameter u satisfying:
    //   Resolvent cubic: u³ - p·u² - 4r·u + (4pr - q²) = 0
    double cubic_a = -p;
    double cubic_b = -4.0 * r;
    double cubic_c = 4.0 * p * r - q * q;

    double cubic_roots[3];
    int num_cubic = solve_cubic(cubic_a, cubic_b, cubic_c, cubic_roots);
    printf("DEBUG: Cubic solver found %d roots\n", num_cubic);
    for (int i = 0; i < num_cubic; i++) {
        printf("  cubic_root[%d] = %.6f\n", i, cubic_roots[i]);
    }

    // Select resolvent root with largest u² = s - p for numerical stability
    //
    // CRITICAL: The constraint is s > p (NOT s > 0) for u² to be positive.
    // When p is negative (common case), negative values of s can be valid.
    // We choose the root that maximizes u² for best numerical conditioning.
    //
    // Example: If p = -1000 and s = -500, then u² = -500 - (-1000) = 500 > 0 ✓
    double s = cubic_roots[0];
    double max_u_squared = s - p;

    for (int i = 1; i < num_cubic; i++) {
        double u_squared_candidate = cubic_roots[i] - p;
        if (u_squared_candidate > max_u_squared) {
            s = cubic_roots[i];
            max_u_squared = u_squared_candidate;
        }
    }

    // Verify we have a valid root
    if (max_u_squared <= 0) {
        // This should never happen for a valid quartic with real roots
        fprintf(
            stderr, "Error: No valid resolvent root found (all give u² ≤ 0)\n"
        );
        return 0;
    }

    // Factor quartic into two quadratics
    double s2 = s * s;
    double u_squared = s - p;

    if (u_squared < 0) {
        // Numerical issue - fallback or no real roots
        return 0;
    }

    double u = sqrt(u_squared);

    // Factor depressed quartic using P² - Q² form:
    //   (y² + u/2)² - [(u-p)y² - qy + (u²/4 - r)] = 0
    //
    // This factors as (P + Q)(P - Q) = 0 where:
    //   P = y² + u/2
    //   Q = u·y - q/(2u)  with u = √(u² = s - p)
    //
    // Expanding P ± Q gives two quadratics:
    //   Q1: y² + u·y + (s/2 - q/(2u)) = 0
    //   Q2: y² - u·y + (s/2 + q/(2u)) = 0

    // Handle case where u ≈ 0
    if (fabs(u) < 1e-12) {
        // Both quadratics are identical: y² + s = 0
        if (s < 0) {
            double y = sqrt(-s);
            roots[0] = y - a / 4.0;
            roots[1] = -y - a / 4.0;
            return 2;
        }
        return 0;
    }

    double e1 = s / 2.0 - q / (2.0 * u);
    double e2 = s / 2.0 + q / (2.0 * u);
    printf("DEBUG: s² = %.6f, u² = %.6f\n", s2, u_squared);
    printf("DEBUG: Quadratic Q1: y² + %.6f·y + %.6f = 0\n", u, e1);
    printf("DEBUG: Quadratic Q2: y² - %.6f·y + %.6f = 0\n", u, e2);

    // Solve first quadratic: y² + u*y + e1 = 0
    int count = 0;
    double disc1 = u * u - 4.0 * e1;
    if (disc1 >= 0) {
        double sqrt_disc1 = sqrt(disc1);
        double y1 = (-u + sqrt_disc1) / 2.0;
        double y2 = (-u - sqrt_disc1) / 2.0;
        printf(
            "DEBUG: Q1 discriminant: %.6f, y-roots: %.6f, %.6f\n", disc1, y1, y2
        );
        roots[count++] = (-u + sqrt_disc1) / 2.0 - a / 4.0;
        roots[count++] = (-u - sqrt_disc1) / 2.0 - a / 4.0;
    }

    // Solve second quadratic: y² - u*y + e2 = 0
    double disc2 = u * u - 4.0 * e2;
    if (disc2 >= 0) {
        double sqrt_disc2 = sqrt(disc2);
        double y3 = (u + sqrt_disc2) / 2.0;
        double y4 = (u - sqrt_disc2) / 2.0;
        printf(
            "DEBUG: Q2 discriminant: %.6f, y-roots: %.6f, %.6f\n", disc2, y3, y4
        );
        roots[count++] = (u + sqrt_disc2) / 2.0 - a / 4.0;
        roots[count++] = (u - sqrt_disc2) / 2.0 - a / 4.0;
    }

    printf("DEBUG: Transform offset: a/4 = %.6f\n", a / 4.0);
    printf("DEBUG: Final λ roots: ");
    for (int i = 0; i < count; i++) {
        printf("%.6f ", roots[i]);
    }
    printf("\n");

    printf("DEBUG: Quartic roots before sort:\n");
    for (int i = 0; i < count; i++) {
        printf("  root[%d] = %.6f\n", i, roots[i]);
    }

    // Sort roots in descending order (largest first)
    for (int i = 0; i < count - 1; i++) {
        for (int j = i + 1; j < count; j++) {
            if (roots[j] > roots[i]) {
                double temp = roots[i];
                roots[i] = roots[j];
                roots[j] = temp;
            }
        }
    }

    return count;
}

/** Compute eigenvector for given eigenvalue using cofactor method
 *
 * For symmetric matrix N with eigenvalue λ, we solve (N - λI)v = 0.
 * The eigenvector v is any non-zero column of the cofactor matrix adj(N - λI).
 *
 * Algorithm:
 *   1. Form M = N - λI (will have rank 3, one zero eigenvalue)
 *   2. Compute all 4 columns of the cofactor matrix of M
 *   3. Select the column with largest norm (for numerical stability)
 *   4. Normalize to unit length
 *
 * NOTE: Mat4d column-major structure means each cofactor column is computed
 * from the appropriate det3() calls using the correct row/column elements.
 *
 * @returns Normalized eigenvector as a quaternion [x, y, z, w]
 */
static Quatd compute_eigenvector_for_eigenvalue(Mat4d N, double lambda) {
    // Form N - λI
    Mat4d M = N;
    M.xi -= lambda;
    M.yj -= lambda;
    M.zk -= lambda;
    M.tw -= lambda;

    // Since rank(M) should be 3 (one zero eigenvalue),
    // one row is linearly dependent on others
    // Find null space by solving M·v = 0

    Quatd v = {0, 0, 0, 0};

    // Use the largest column from cofactor matrix
    // The eigenvector is any non-zero column of adj(M) where adj = cofactor
    // matrix

    // For numerical stability, find the column with largest norm

    // Column 0 cofactors
    double c0x = det3(M.yj, M.zj, M.tj, M.yk, M.zk, M.tk, M.yw, M.zw, M.tw);
    double c0y = -det3(M.yi, M.zi, M.ti, M.yk, M.zk, M.tk, M.yw, M.zw, M.tw);
    double c0z = det3(M.yi, M.zi, M.ti, M.yj, M.zj, M.tj, M.yw, M.zw, M.tw);
    double c0w = -det3(M.yi, M.zi, M.ti, M.yj, M.zj, M.tj, M.yk, M.zk, M.tk);

    Vec4d col0 = {c0x, c0y, c0z, c0w};
    double norm0_sq = c0x * c0x + c0y * c0y + c0z * c0z + c0w * c0w;

    // Column 1 cofactors
    double c1x = -det3(M.xj, M.zj, M.tj, M.xk, M.zk, M.tk, M.xw, M.zw, M.tw);
    double c1y = det3(M.xi, M.zi, M.ti, M.xk, M.zk, M.tk, M.xw, M.zw, M.tw);
    double c1z = -det3(M.xi, M.zi, M.ti, M.xj, M.zj, M.tj, M.xw, M.zw, M.tw);
    double c1w = det3(M.xi, M.zi, M.ti, M.xj, M.zj, M.tj, M.xk, M.zk, M.tk);

    Vec4d col1 = {c1x, c1y, c1z, c1w};
    double norm1_sq = c1x * c1x + c1y * c1y + c1z * c1z + c1w * c1w;

    // Column 2 cofactors
    double c2x = det3(M.xj, M.yj, M.tj, M.xk, M.yk, M.tk, M.xw, M.yw, M.tw);
    double c2y = -det3(M.xi, M.yi, M.ti, M.xk, M.yk, M.tk, M.xw, M.yw, M.tw);
    double c2z = det3(M.xi, M.yi, M.ti, M.xj, M.yj, M.tj, M.xw, M.yw, M.tw);
    double c2w = -det3(M.xi, M.yi, M.ti, M.xj, M.yj, M.tj, M.xk, M.yk, M.tk);

    Vec4d col2 = {c2x, c2y, c2z, c2w};
    double norm2_sq = c2x * c2x + c2y * c2y + c2z * c2z + c2w * c2w;

    // Column 3 cofactors
    double c3x = -det3(M.xj, M.yj, M.zj, M.xk, M.yk, M.zk, M.xw, M.yw, M.zw);
    double c3y = det3(M.xi, M.yi, M.zi, M.xk, M.yk, M.zk, M.xw, M.yw, M.zw);
    double c3z = -det3(M.xi, M.yi, M.zi, M.xj, M.yj, M.zj, M.xw, M.yw, M.zw);
    double c3w = det3(M.xi, M.yi, M.zi, M.xj, M.yj, M.zj, M.xk, M.yk, M.zk);

    Vec4d col3 = {c3x, c3y, c3z, c3w};
    double norm3_sq = c3x * c3x + c3y * c3y + c3z * c3z + c3w * c3w;

    // Pick column with largest norm (most numerically stable)
    v = col0;
    double max_norm_sq = norm0_sq;

    if (norm1_sq > max_norm_sq) {
        v = col1;
        max_norm_sq = norm1_sq;
    }
    if (norm2_sq > max_norm_sq) {
        v = col2;
        max_norm_sq = norm2_sq;
    }
    if (norm3_sq > max_norm_sq) {
        v = col3;
        max_norm_sq = norm3_sq;
    }

    // Normalize
    return Quatd_normalize(v);
}

/** Ferrari's Method for solving quartic characteristic polynomial
 *
 * Reference: https://mathworld.wolfram.com/QuarticEquation.html
 *
 * Steps:
 *   1. Compute characteristic polynomial det(N - λI) = 0 via determinants
 *   2. Depress the quartic (eliminate λ³ term) by substituting λ = y - a/4
 *   3. Solve MathWorld resolvent cubic: u³ - p·u² - 4r·u + (4pr - q²) = 0
 *   4. Select resolvent root u with largest u² = u - p (ensures u² > 0)
 *   5. Factor depressed quartic into two quadratics using u
 *   6. Solve both quadratics to find all 4 eigenvalues
 *   7. Select MOST POSITIVE eigenvalue (per Horn 1987 paper)
 *   8. Compute corresponding eigenvector using cofactor method
 *
 * @returns RETURN_OK on success, or falls back to power iteration on failure
 */
static ReturnCode ferrari(const Mat4d N, Quatd* q_vec) {
    double poly_coeffs[5] = {0};
    // Fourth-order term
    poly_coeffs[0] = 1.0;
    // Third-order term: -trace(N)
    poly_coeffs[1] = -(N.xi + N.yj + N.zk + N.tw);

    // Second-order term: sum of 2×2 principal minors
    // These are determinants of 2×2 submatrices from diagonal element pairs
    //
    // NOTE: Due to column-major Mat4d, diagonal elements are: xi, yj, zk, tw
    double m01 = det2(N.xi, N.yi, N.xj, N.yj);
    double m02 = det2(N.xi, N.zi, N.xk, N.zk);
    double m03 = det2(N.xi, N.ti, N.xw, N.tw);
    double m12 = det2(N.yj, N.zj, N.yk, N.zk);
    double m13 = det2(N.yj, N.tj, N.yw, N.tw);
    double m23 = det2(N.zk, N.tk, N.zw, N.tw);
    poly_coeffs[2] = m01 + m02 + m03 + m12 + m13 + m23;
    printf(
        "2x2 minors: m01=%.6f, m02=%.6f, m03=%.6f, m12=%.6f, m13=%.6f, "
        "m23=%.6f\n",
        m01,
        m02,
        m03,
        m12,
        m13,
        m23
    );
    printf("Sum of 2x2 minors: %.6f\n", m01 + m02 + m03 + m12 + m13 + m23);

    // First-order term: -sum of 3×3 principal minors
    // These are determinants when one row/col is removed
    double n0 = det3(N.yj, N.zj, N.tj, N.yk, N.zk, N.tk, N.yw, N.zw, N.tw);
    double n1 = det3(N.xi, N.zi, N.ti, N.xk, N.zk, N.tk, N.xw, N.zw, N.tw);
    double n2 = det3(N.xi, N.yi, N.ti, N.xj, N.yj, N.tj, N.xw, N.yw, N.tw);
    double n3 = det3(N.xi, N.yi, N.zi, N.xj, N.yj, N.zj, N.xk, N.yk, N.zk);
    poly_coeffs[3] = -(n0 + n1 + n2 + n3);
    printf("3x3 minors: n0=%.6f, n1=%.6f, n2=%.6f, n3=%.6f\n", n0, n1, n2, n3);
    printf("Sum of 3x3 minors: %.6f\n", n0 + n1 + n2 + n3);

    // Coefficient of λ⁰: det(N) - full 4×4 determinant
    // NOTE: For symmetric matrices (like Horn's N), row/column expansion are
    // equivalent. We use column 0 expansion because it's natural with
    // column-major storage.
    double M_01 = det3(N.yi, N.yk, N.yw, N.zi, N.zk, N.zw, N.ti, N.tk, N.tw);
    double M_02 = det3(N.yi, N.yj, N.yw, N.zi, N.zj, N.zw, N.ti, N.tj, N.tw);
    double M_03 = det3(N.yi, N.yj, N.yk, N.zi, N.zj, N.zk, N.ti, N.tj, N.tk);
    poly_coeffs[4] = N.xi * n0 - N.xj * M_01 + N.xk * M_02 - N.xw * M_03;
    printf(
        "Determinant minors: M_01=%.6f, M_02=%.6f, M_03=%.6f\n",
        M_01,
        M_02,
        M_03
    );
    printf(
        "det(N) = %.6f * %.6f - %.6f * %.6f + %.6f * %.6f - %.6f * %.6f = "
        "%.6f\n",
        N.xi,
        n0,
        N.xj,
        M_01,
        N.xk,
        M_02,
        N.xw,
        M_03,
        poly_coeffs[4]
    );

    printf(
        "Characteristic polynomial: %.6fλ⁴ + %.6fλ³ + %.6fλ² + %.6fλ + %.6f\n",
        poly_coeffs[0],
        poly_coeffs[1],
        poly_coeffs[2],
        poly_coeffs[3],
        poly_coeffs[4]
    );

    double eigenvalues[4];
    uint32_t num_roots = solve_quartic(
        poly_coeffs[1],
        poly_coeffs[2],
        poly_coeffs[3],
        poly_coeffs[4],
        eigenvalues
    );

    printf("Found %d real eigenvalues:", num_roots);
    for (uint32_t i = 0; i < num_roots; i++) {
        printf(" %.6f", eigenvalues[i]);
    }
    printf("\n");

    // Validation check
    if (num_roots == 0) {
        fprintf(
            stderr,
            "Warning: No real eigenvalues found (unexpected for symmetric "
            "matrix)\n"
        );
        fprintf(stderr, "Falling back to power iteration\n");
        // Fall back to power iteration
        return power_iteration(N, MAX_ITERATIONS, q_vec);
    }

    printf("All eigenvalues: ");
    for (uint32_t i = 0; i < num_roots; i++) {
        printf("%.6f ", eigenvalues[i]);
    }
    printf("\n");

    // Select the MOST POSITIVE eigenvalue (per Horn 1987)
    //
    // Horn's method requires the eigenvector corresponding to the most positive
    // eigenvalue (not largest absolute value, not largest magnitude).
    // For point cloud alignment, this gives the optimal rotation quaternion.
    double lambda_max = eigenvalues[0];
    uint32_t max_abs_idx = 0;
    for (uint32_t i = 1; i < num_roots; i++) {
        double lambda_new = eigenvalues[i];
        if (lambda_new > lambda_max) {
            max_abs_idx = i;
            lambda_max = lambda_new;
        }
    }
    printf(
        "Ferrari Eigenvalue: %.6f (index %d)\n", lambda_max, max_abs_idx
    );

    // Degeneracy detection
    if (num_roots >= 2) {
        double ratio = eigenvalues[1] / eigenvalues[0];
        if (ratio > 0.95) {
            printf(
                "Warning: Near-degenerate eigenvalues (ratio=%.6f)\n", ratio
            );
        }
    }

    // Compute eigenvector
    *q_vec = compute_eigenvector_for_eigenvalue(N, lambda_max);
    return RETURN_OK;
}

/* Calculate an alignment transform using Horn's algorithm
 *
 * Reference: https://people.csail.mit.edu/bkph/papers/Absolute_Orientation.pdf
 */
Transform3d align_points(
    const Vec3d* source, const Vec3d* target, uint32_t n_points
) {
    // Handle edge case
    if (n_points == 0) {
        Transform3d t = {
            .translation = VEC3D_ZERO,
            .rotation = QUATD_IDENTITY,
        };
        return t;
    }

    // Step 1: Compute centroids
    Vec3d centroid_source = compute_centroid(source, n_points);
    printf(
        "Source centroid: [%12.6f, %12.6f, %12.6f]\n",
        centroid_source.x,
        centroid_source.y,
        centroid_source.z
    );
    Vec3d centroid_target = compute_centroid(target, n_points);
    printf(
        "Target centroid: [%12.6f, %12.6f, %12.6f]\n",
        centroid_target.x,
        centroid_target.y,
        centroid_target.z
    );

    // Step 2: Build cross-covariance matrix (3x3)
    double Sxx = 0.0, Sxy = 0.0, Sxz = 0.0;
    double Syx = 0.0, Syy = 0.0, Syz = 0.0;
    double Szx = 0.0, Szy = 0.0, Szz = 0.0;
    for (uint32_t i = 0; i < n_points; i++) {
        Vec3d src = Vec3d_sub(source[i], centroid_source);
        Vec3d tgt = Vec3d_sub(target[i], centroid_target);

        Sxx += src.x * tgt.x;
        Sxy += src.x * tgt.y;
        Sxz += src.x * tgt.z;
        Syx += src.y * tgt.x;
        Syy += src.y * tgt.y;
        Syz += src.y * tgt.z;
        Szx += src.z * tgt.x;
        Szy += src.z * tgt.y;
        Szz += src.z * tgt.z;
    }

    /* Step 3: Build Horn's 4×4 symmetric matrix N
     *
     * Horn's 1987 paper uses quaternion order [w, x, y, z] (scalar-first).
     * We use [x, y, z, w] (scalar-last), so we rearrange the matrix
     * accordingly.
     *
     * N matrix layout for [x, y, z, w] quaternion convention:
     *
     *       col 0 (x)      col 1 (y)      col 2 (z)      col 3 (w)
     *   ┌────────────────────────────────────────────────────────────┐
     * 0 │ Sxx-Syy-Szz    Sxy+Syx        Szx+Sxz        Syz-Szy       │
     * 1 │ Sxy+Syx       -Sxx+Syy-Szz    Syz+Szy        Szx-Sxz       │
     * 2 │ Szx+Sxz        Syz+Szy       -Sxx-Syy+Szz    Sxy-Syx       │
     * 3 │ Syz-Szy        Szx-Sxz        Sxy-Syx        Sxx+Syy+Szz   │
     *   └────────────────────────────────────────────────────────────┘
     *
     * This is a symmetric matrix (N = Nᵀ) where the eigenvector corresponding
     * to the MOST POSITIVE eigenvalue gives the optimal rotation quaternion.
     */
    Mat4d N;
    // Column 0 (x component of quaternion)
    N.xi = Sxx - Syy - Szz;
    N.xj = Sxy + Syx;
    N.xk = Szx + Sxz;
    N.xw = Syz - Szy;

    // Column 1 (y component)
    N.yi = Sxy + Syx;
    N.yj = -Sxx + Syy - Szz;
    N.yk = Syz + Szy;
    N.yw = Szx - Sxz;

    // Column 2 (z component)
    N.zi = Szx + Sxz;
    N.zj = Syz + Szy;
    N.zk = -Sxx - Syy + Szz;
    N.zw = Sxy - Syx;

    // Column 3 (w component)
    N.ti = Syz - Szy;
    N.tj = Szx - Sxz;
    N.tk = Sxy - Syx;
    N.tw = Sxx + Syy + Szz;

    printf("N matrix:\n");
    Mat4d_print(N);

    printf("\nTest case N = diag(2, 3, 5, 7)...\n");
    Mat4d N_test = MAT4D_IDENTITY;
    N_test.xi = 2.0;
    N_test.yj = 3.0;
    N_test.zk = 5.0;
    N_test.tw = 7.0;
    printf("Test case Ferrari method...\n");
    Quatd q_test = QUATD_IDENTITY;
    ferrari(N_test, &q_test);

    printf("\nContinuing with provided N...\n");
    // Find largest diagonal element
    double diag[4] = {N.xi, N.yj, N.zk, N.tw};
    int max_idx = 0;
    double max_val = diag[0];
    for (int i = 1; i < 4; i++) {
        if (diag[i] > max_val) {
            max_val = diag[i];
            max_idx = i;
        }
    }

    // Step 4: Eigenvalue solver

    Vec4d q_vec;
    // q_vec = (Quatd)QUATD_IDENTITY;
    switch (max_idx) {
        case 0:
            q_vec = (Vec4d){N.xi, N.xj, N.xk, N.xw};
            break;
        case 1:
            q_vec = (Vec4d){N.yi, N.yj, N.yk, N.yw};
            break;
        case 2:
            q_vec = (Vec4d){N.zi, N.zj, N.zk, N.zw};
            break;
        case 3:
            q_vec = (Vec4d){N.ti, N.tj, N.tk, N.tw};
            break;
    }
    q_vec = Quatd_normalize(q_vec);

#if USE_POWER_ITERATION
    printf("\nPower Iteration method...\n");
    ReturnCode rc_power = power_iteration(N, MAX_ITERATIONS, &q_vec);
    assert(rc_power != RETURN_ITERATION_FAIL && "Convergence NOT reached!\n");
#else
    printf("\nFerrari method...\n");
    ReturnCode rc_ferrari = ferrari(N, &q_vec);
    assert(rc_ferrari != RETURN_FAIL && "Convergence NOT reached!\n");
#endif

    printf("\nBelow results are for Ferrari method...\n");

    // Ensure the quaternion is normalized
    double q_mag = Quatd_length(q_vec);
    if (fabs(q_mag - 1.0) > 1.0e-9) {
        q_vec = Quatd_normalize(q_vec);
    }

    // Ensure consistent quaternion convention
    if (q_vec.w < 0.0) {
        q_vec.x = -q_vec.x;
        q_vec.y = -q_vec.y;
        q_vec.z = -q_vec.z;
        q_vec.w = -q_vec.w;
    }

    // Step 5: Compute translation
    Vec3d rotated_source_centroid = Quatd_rotate_vec3d(q_vec, centroid_source);
    Vec3d translation = Vec3d_sub(centroid_target, rotated_source_centroid);
    Transform3d t = {
        .translation = translation,
        .rotation = q_vec,
    };
    return t;
}

Transform3d align_points_from_maps(Vec3dMap* measured, Vec3dMap* nominals) {
    Transform3d default_xform = {
        .translation = VEC3D_ZERO,
        .rotation = QUATD_IDENTITY,
    };

    // Nothing to do
    if (measured->count == 0 || nominals->count == 0) {
        return default_xform;
    }

    // Allocate worst-case size (all measured points might match)
    Vec3d* measured_array = malloc(measured->count * sizeof(Vec3d));
    Vec3d* nominal_array = malloc(measured->count * sizeof(Vec3d));

    if (!measured_array || !nominal_array) {
        // Handle allocation failure
        free(measured_array);
        free(nominal_array);
        return default_xform;
    }

    // Iterate through measured points and find corresponding nominals
    uint32_t n_matched = 0;
    for (uint32_t i = 0; i < measured->capacity; i++) {
        Vec3dEntry* entry = measured->entries[i];
        while (entry) {
            // Look up this key in nominals map
            Vec3d* nominal_pt = Vec3dMap_get(nominals, entry->key);

            if (nominal_pt != NULL) {
                // Found a match - add to alignment arrays
                measured_array[n_matched] = entry->value;
                nominal_array[n_matched] = *nominal_pt;
                n_matched++;
            }

            // Handle hash collisions
            entry = entry->next;
        }
    }

    Transform3d result;
    if (n_matched >= 3) {
        // Need at least 3 points for 3D alignment
        result = align_points(measured_array, nominal_array, n_matched);
    } else {
        // Not enough matches
        result = default_xform;
    }

    free(measured_array);
    free(nominal_array);

    return result;
}
