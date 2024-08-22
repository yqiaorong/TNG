#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>

void create_mask(float *Coordinates, bool *mask, int num_coords, float xmin, float xmax, float ymin, float ymax, float zmin, float zmax) {
    #pragma omp parallel for
    for (int i = 0; i < num_coords; i++) {
        float x = Coordinates[3*i];
        float y = Coordinates[3*i+1];
        float z = Coordinates[3*i+2];
        if (xmin <= x && x <= xmax && ymin <= y && y <= ymax && zmin <= z && z <= zmax) {
            mask[i] = true;
        } else {
            mask[i] = false;
        }
    }
}

int select_sub_coords(float *Coordinates, float *selected_coords, bool *mask, int num_coords) {
    int idx = 0;
    for (int i = 0; i < num_coords; i++) {
        if (mask[i]) {
            selected_coords[3*idx] = Coordinates[3*i];
            selected_coords[3*idx+1] = Coordinates[3*i+1];
            selected_coords[3*idx+2] = Coordinates[3*i+2];
            idx++;
        }
    }
    return idx;
}