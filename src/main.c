#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "point_align.h"

void test_vec3_hashmap(void) {
    printf("=== Testing Vec3Map ===\n");

    // Initialize map
    Vec3dMap map;
    Vec3dMap_init(&map);

    // Add some points
    Vec3dMap_put(&map, "point_a", Vec3d_init(1.0f, 2.0f, 3.0f));
    Vec3dMap_put(&map, "point_b", Vec3d_init(4.0f, 5.0f, 6.0f));
    Vec3dMap_put(&map, "point_c", Vec3d_init(7.0f, 8.0f, 9.0f));

    printf("Added %u points\n", map.count);

    // Retrieve points
    Vec3d *pa = Vec3dMap_get(&map, "point_a");
    Vec3d *pb = Vec3dMap_get(&map, "point_b");
    Vec3d *pc = Vec3dMap_get(&map, "point_c");
    Vec3d *missing = Vec3dMap_get(&map, "point_d");

    if (pa) printf("point_a: {%.1f, %.1f, %.1f}\n", pa->x, pa->y, pa->z);
    if (pb) printf("point_b: {%.1f, %.1f, %.1f}\n", pb->x, pb->y, pb->z);
    if (pc) printf("point_c: {%.1f, %.1f, %.1f}\n", pc->x, pc->y, pc->z);
    printf(
        "point_d (missing): %s\n", missing ? "FOUND (ERROR)" : "NULL (correct)"
    );

    // Update a value
    Vec3dMap_put(&map, "point_a", Vec3d_init(10.0f, 11.0f, 12.0f));
    pa = Vec3dMap_get(&map, "point_a");
    printf("Updated point_a: {%.1f, %.1f, %.1f}\n", pa->x, pa->y, pa->z);
    printf("Count after update: %u (should still be 3)\n", map.count);

    // Clean up
    Vec3dMap_free(&map);
    printf("Freed map\n\n");
}

double rand_range_f64(float min, float max) {
    double result = min + ((float)rand() / (float)RAND_MAX) * (max - min);
    return result;
}

int32_t rand_range_i32(int32_t min, int32_t max) {
    int32_t result = min + (rand() % (max - min + 1));
    return result;
}

void test_point_alignment(void) {
    static const uint32_t POINT_COUNT = 5;
    Vec3d nominal[POINT_COUNT];
    for (uint32_t i = 0; i < POINT_COUNT; ++i) {
        Vec3d *nom = &nominal[i];
        nom->x = rand_range_f64(-100.0, 100.0);
        nom->y = rand_range_f64(-100.0, 100.0);
        nom->z = rand_range_f64(-100.0, 100.0);
    }

    Transform3d xform = {
        .translation =
            {
                .x = rand_range_f64(-100.0, 100.0),
                .y = rand_range_f64(-100.0, 100.0),
                .z = rand_range_f64(-100.0, 100.0),
            },
        .rotation = Quatd_normalize((Quatd){
            .x = rand_range_f64(-1.0, 1.0),
            .y = rand_range_f64(-1.0, 1.0),
            .z = rand_range_f64(-1.0, 1.0),
            .w = rand_range_f64(-1.0, 1.0),
        }),
    };
    // Ensure consistent quaternion convention
    if (xform.rotation.w < 0.0) {
        xform.rotation.x = -xform.rotation.x;
        xform.rotation.y = -xform.rotation.y;
        xform.rotation.z = -xform.rotation.z;
        xform.rotation.w = -xform.rotation.w;
    }
    printf("Transformation:\n");
    printf(
        "  Translation:    [%12.6f, %12.6f, %12.6f]\n",
        xform.translation.x,
        xform.translation.y,
        xform.translation.z
    );
    printf(
        "  Rotation:       [%12.6f, %12.6f, %12.6f, %12.6f]\n",
        xform.rotation.x,
        xform.rotation.y,
        xform.rotation.z,
        xform.rotation.w
    );

    Vec3d measured[POINT_COUNT];
    for (uint32_t i = 0; i < POINT_COUNT; ++i) {
        Vec3d pt = transform_point(nominal[i], xform);
        Vec3d *meas = &measured[i];
        meas->x = pt.x;
        meas->y = pt.y;
        meas->z = pt.z;
    }

    Transform3d alignment = align_points(nominal, measured, POINT_COUNT);

    printf("\nComparison:\n");
    printf(
        "  Applied Translation:    [%12.6f, %12.6f, %12.6f]\n",
        xform.translation.x,
        xform.translation.y,
        xform.translation.z
    );
    printf(
        "  Calculated Translation: [%12.6f, %12.6f, %12.6f]\n",
        alignment.translation.x,
        alignment.translation.y,
        alignment.translation.z
    );
    printf(
        "  Deviation:              [%12.6f, %12.6f, %12.6f]\n",
        alignment.translation.x - xform.translation.x,
        alignment.translation.y - xform.translation.y,
        alignment.translation.z - xform.translation.z
    );
    printf(
        "\n  Applied Rotation:       [%12.6f, %12.6f, %12.6f, %12.6f]\n",
        xform.rotation.x,
        xform.rotation.y,
        xform.rotation.z,
        xform.rotation.w
    );
    printf(
        "  Calculated Rotation:    [%12.6f, %12.6f, %12.6f, %12.6f]\n",
        alignment.rotation.x,
        alignment.rotation.y,
        alignment.rotation.z,
        alignment.rotation.w
    );
    printf(
        "  Deviation:              [%12.6f, %12.6f, %12.6f, %12.6f]\n",
        alignment.rotation.x - xform.rotation.x,
        alignment.rotation.y - xform.rotation.y,
        alignment.rotation.z - xform.rotation.z,
        alignment.rotation.w - xform.rotation.w
    );

    double rmsd = compute_rmsd(nominal, measured, POINT_COUNT, alignment);
    printf("\nRMSD: %12.6f\n", rmsd);
    assert(rmsd < 1e-9);
}

int main(void) {
    test_vec3_hashmap();
    srand(time(NULL));
    for (uint32_t i = 0; i < 1000; ++i) {
        printf(
            "---------------------------------------------------------------\n"
        );
        printf("TEST %d\n", i);
        printf(
            "---------------------------------------------------------------\n"
        );
        test_point_alignment();
    }
    return 0;
}
