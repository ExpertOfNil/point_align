#include <assert.h>
#include <stdio.h>

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

int main(void) {
    test_vec3_hashmap();
    return 0;
}
