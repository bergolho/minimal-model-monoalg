//
// Created by sachetto on 16/03/2021.
//
#include "domain_helpers.h"
#include "custom_mesh_info_data.h"

#include "../3dparty/sds/sds.h"
#include "../3dparty/stb_ds.h"
#include "../config_helpers/config_helpers.h"
#include "../logger/logger.h"
#include "../utils/utils.h"
#include <time.h>
#include <unistd.h>
#include <float.h>

static void set_human_mesh_fibrosis_from_file(struct grid *grid, int type, const char *filename, int size) {

    FILE *file = fopen(filename, "r");

    if(!file) {
        printf("Error opening file %s!!\n", filename);
        exit(0);
    }

    real_cpu **scar_mesh = (real_cpu **)malloc(sizeof(real_cpu *) * size);
    for(int i = 0; i < size; i++) {
        scar_mesh[i] = (real_cpu *)malloc(sizeof(real_cpu) * 3);
        if(scar_mesh[i] == NULL) {
            printf("Failed to allocate memory\n");
            exit(0);
        }
    }
    real_cpu dummy1, dummy2; // unused values

    int i = 0;

    while(!feof(file)) {
        fscanf(file, "%lf,%lf,%lf,%lf,%lf\n", &scar_mesh[i][0], &scar_mesh[i][1], &scar_mesh[i][2], &dummy1, &dummy2);
        i++;
    }

    fclose(file);

    sort_vector(scar_mesh, size);

    struct cell_node *grid_cell = grid->first_cell;
    while(grid_cell != 0) {

        real_cpu center_x = grid_cell->center.x;
        real_cpu center_y = grid_cell->center.y;
        real_cpu center_z = grid_cell->center.z;

        if((grid_cell->discretization.x == 100.0) && (DHZB_MESH_TISSUE_TYPE(grid_cell) == type)) {
            int index = inside_mesh(scar_mesh, center_x, center_y, center_z, 0, size - 1);
            grid_cell->active = (index != -1);
        }

        grid_cell = grid_cell->next;
    }

    for(int k = 0; k < size; k++) {
        free(scar_mesh[k]);
    }

    free(scar_mesh);
}

static int set_human_mesh(struct grid *the_grid, struct config *config, int n_custom_data, set_custom_data_for_mesh_fn set_custom_data) {
    real_cpu original_discretization = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, original_discretization, config, "original_discretization");

    real_cpu start_discretization = original_discretization;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, start_discretization, config, "start_discretization");

    real_cpu max_discretization = start_discretization;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, max_discretization, config, "max_discretization");

    uint32_t size = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(uint32_t, size, config, "num_volumes");

    char *mesh_file;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(mesh_file, config, "mesh_file");


    if(original_discretization != 800 || original_discretization != 500) {
        log_error_and_exit("Invalid original_discretization: %lf. Valid values are 500 or 800\n", original_discretization);
    }

    uint64_t num_loaded = set_custom_mesh_from_file(the_grid, mesh_file, size, original_discretization, n_custom_data, set_custom_data);

    log_info("Read %d volumes from file: %s\n", num_loaded, mesh_file);

    free(mesh_file);

    the_grid->start_discretization = SAME_POINT3D(start_discretization);
    the_grid->max_discretization = SAME_POINT3D(max_discretization);

    int remaining_refinements = (int)((original_discretization / start_discretization) - 1);
    refine_grid(the_grid, remaining_refinements);

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_human_mesh) {
    return set_human_mesh(the_grid, config, 0, NULL);
}

static void refine_fibrotic_cells(struct grid *the_grid) {

    assert(the_grid);

    struct cell_node *grid_cell, *auxiliar_grid_cell;

    struct dhzb_mesh_info *mesh_info;

    grid_cell = the_grid->first_cell;
    while(grid_cell != 0) {

        mesh_info = DHZB_MESH_INFO(grid_cell);

        if(grid_cell->active && mesh_info->tissue_type == SCAR) {
            auxiliar_grid_cell = grid_cell;
            grid_cell = grid_cell->next;
            refine_cell(auxiliar_grid_cell, NULL, NULL);
            the_grid->number_of_cells += 7;
        } else {
            grid_cell = grid_cell->next;
        }
    }
}

static void refine_border_zone_cells(struct grid *the_grid) {

    assert(the_grid);

    struct cell_node *grid_cell, *auxiliar_grid_cell;

    struct dhzb_mesh_info *mesh_info;

    grid_cell = the_grid->first_cell;
    while(grid_cell != 0) {

        mesh_info = DHZB_MESH_INFO(grid_cell);

        if(grid_cell->active && mesh_info->tissue_type == BZ) {
            auxiliar_grid_cell = grid_cell;
            grid_cell = grid_cell->next;
            refine_cell(auxiliar_grid_cell, NULL, NULL);
            the_grid->number_of_cells += 7;
        } else {
            grid_cell = grid_cell->next;
        }
    }
}

static void set_human_mesh_fibrosis(struct grid *grid, real_cpu phi, unsigned seed, real_cpu big_scar_center_x,
                             real_cpu big_scar_center_y, real_cpu big_scar_center_z, real_cpu small_scar_center_x,
                             real_cpu small_scar_center_y, real_cpu small_scar_center_z) {

    if(seed == 0)
        seed = (unsigned)time(NULL) + getpid();

    srand(seed);

    log_info("Using %u as seed\n", seed);

    real_cpu bz_size_big = 0;
    real_cpu bz_size_small = 0;
    real_cpu dist_big = 0;
    real_cpu dist_small = 0;

    log_info("Calculating fibrosis using phi: %lf\n", phi);
    struct cell_node *grid_cell = grid->first_cell;

    while(grid_cell != NULL) {

        if(grid_cell->active) {
            if(DHZB_MESH_TISSUE_TYPE(grid_cell) == SCAR) {
                grid_cell->can_change = false;
                real_cpu p = (real_cpu)(rand()) / (RAND_MAX);
                if(p < phi)
                    grid_cell->active = false;
            } else if(DHZB_MESH_TISSUE_TYPE(grid_cell) == BZ) {
                real_cpu centerX = grid_cell->center.x;
                real_cpu centerY = grid_cell->center.y;
                real_cpu centerZ = grid_cell->center.z;
                if(DHZB_MESH_LOCATION(grid_cell) == BIG_SCAR) {
                    dist_big = sqrt((centerX - big_scar_center_x) * (centerX - big_scar_center_x) +
                                    (centerY - big_scar_center_y) * (centerY - big_scar_center_y) +
                                    (centerZ - big_scar_center_z) * (centerZ - big_scar_center_z));
                    if(dist_big > bz_size_big) {
                        bz_size_big = dist_big;
                    }
                } else if(DHZB_MESH_LOCATION(grid_cell) == SMALL_SCAR) {
                    dist_small = sqrt((centerX - small_scar_center_x) * (centerX - small_scar_center_x) +
                                      (centerY - small_scar_center_y) * (centerY - small_scar_center_y) +
                                      (centerZ - small_scar_center_z) * (centerZ - small_scar_center_z));
                    if(dist_small > bz_size_small) {
                        bz_size_small = dist_small;
                    }
                }
            }
        }
        grid_cell = grid_cell->next;
    }

    grid_cell = grid->first_cell;
    while(grid_cell != NULL) {

        if(grid_cell->active) {
            if(DHZB_MESH_TISSUE_TYPE(grid_cell) == BZ) {
                real_cpu centerX = grid_cell->center.x;
                real_cpu centerY = grid_cell->center.y;
                real_cpu centerZ = grid_cell->center.z;
                if(DHZB_MESH_LOCATION(grid_cell) == BIG_SCAR) {
                    dist_big = sqrt((centerX - big_scar_center_x) * (centerX - big_scar_center_x) +
                                    (centerY - big_scar_center_y) * (centerY - big_scar_center_y) +
                                    (centerZ - big_scar_center_z) * (centerZ - big_scar_center_z));
                    dist_big = dist_big / bz_size_big;
                    real_cpu phi_local = phi - phi * dist_big;

                    real_cpu p = (real_cpu)(rand()) / (RAND_MAX);
                    if(p < phi_local) {
                        grid_cell->active = false;
                    }
                    grid_cell->can_change = false;
                } else if(DHZB_MESH_LOCATION(grid_cell) == SMALL_SCAR) {
                    dist_small = sqrt((centerX - small_scar_center_x) * (centerX - small_scar_center_x) +
                                      (centerY - small_scar_center_y) * (centerY - small_scar_center_y) +
                                      (centerZ - small_scar_center_z) * (centerZ - small_scar_center_z));
                    dist_small = dist_small / bz_size_small;
                    real_cpu phi_local = phi - phi * dist_small;

                    real_cpu p = (real_cpu)(rand()) / (RAND_MAX);
                    if(p < phi_local) {
                        grid_cell->active = false;
                    }
                    grid_cell->can_change = false;
                }
            }
        }
        grid_cell = grid_cell->next;
    }
}

SET_CUSTOM_DATA_FOR_MESH(set_custom_data_for_dhzb_mesh) {

    INITIALIZE_DHZB_MESH_INFO(cell);

    int id = custom_data[0];

    if(cell->center.x > 79300) {
        DHZB_MESH_LOCATION(cell) = BIG_SCAR;
    } else {
        DHZB_MESH_LOCATION(cell) = SMALL_SCAR;
    }

    if (id == 400) {
        DHZB_MESH_TISSUE_TYPE(cell) = SCAR;
    }
    else if (id == 130 || id == 230 || id == 155 || id == 255) {
        DHZB_MESH_TISSUE_TYPE(cell) = BZ;
    }
    else {
        DHZB_MESH_TISSUE_TYPE(cell) = HEALTH;
        DHZB_MESH_LOCATION(cell) = HEALTH_AREA;
    }

    cell->sigma.fibers.f[0] = custom_data[1];
    cell->sigma.fibers.f[1] = custom_data[2];
    cell->sigma.fibers.f[2] = custom_data[3];

    cell->sigma.fibers.s[0] = custom_data[4];
    cell->sigma.fibers.s[1] = custom_data[5];
    cell->sigma.fibers.s[2] = custom_data[6];

}

SET_SPATIAL_DOMAIN(initialize_grid_with_human_mesh_with_two_scars) {

    int success = set_human_mesh(the_grid, config, 7, set_custom_data_for_dhzb_mesh);

    if(!success) return 0;

    bool fibrotic = false;
    GET_PARAMETER_BOOLEAN_VALUE_OR_USE_DEFAULT(fibrotic, config, "fibrotic");

    if(fibrotic) {

        // Here we refine the scar cells
        refine_fibrotic_cells(the_grid); // 400 um
        refine_fibrotic_cells(the_grid); // 200 um
        refine_fibrotic_cells(the_grid); // 100 um

        // and the border zone
        refine_border_zone_cells(the_grid);
        refine_border_zone_cells(the_grid);
        refine_border_zone_cells(the_grid);

        char *scar_file_big = NULL;
        GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(scar_file_big, config, "big_scar_file");

        char *scar_file_small = NULL;
        GET_PARAMETER_STRING_VALUE_OR_USE_DEFAULT(scar_file_small, config, "small_scar_file");

        if(scar_file_big) {
            log_info("Loading fibrosis patterns from file %s\n", scar_file_big);
            set_human_mesh_fibrosis_from_file(the_grid, 'b', scar_file_big, 2172089);
        }

        if(scar_file_small) {
            log_info("Loading fibrosis patterns from file %s\n", scar_file_small);
            set_human_mesh_fibrosis_from_file(the_grid, 's', scar_file_small, 845051);
        }

        if(!(scar_file_big || scar_file_small)) {

            real_cpu small_scar_center_x = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, small_scar_center_x, config, "small_scar_center_x");

            real_cpu small_scar_center_y = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, small_scar_center_y, config, "small_scar_center_y");

            real_cpu small_scar_center_z = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, small_scar_center_z, config, "small_scar_center_z");

            real_cpu big_scar_center_x = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, big_scar_center_x, config, "big_scar_center_x");

            real_cpu big_scar_center_y = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, big_scar_center_y, config, "big_scar_center_y");

            real_cpu big_scar_center_z = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real, big_scar_center_z, config, "big_scar_center_z");

            real_cpu phi = 0.0;
            GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config, "phi");

            unsigned seed = 0;
            GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, seed, config, "seed");

            log_info("Setting random fibrosis pattern\n");
            set_human_mesh_fibrosis(the_grid, phi, seed, big_scar_center_x, big_scar_center_y, big_scar_center_z, small_scar_center_x, small_scar_center_y,
                                    small_scar_center_z);
        }
    }

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_scar_wedge) {
    char *mesh_file;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(mesh_file, config, "mesh_file");

    char *scar_size;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(scar_size, config, "scar_size");

    real_cpu phi = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, phi, config, "phi");

    unsigned fib_seed = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(unsigned, fib_seed, config, "seed");

    if(fib_seed == 0)
    fib_seed = (unsigned)time(NULL) + getpid();

    char tmp[256];
    sprintf(tmp, "%u", fib_seed);
    shput_dup_value(config->config_data, "seed", "tmp");

    srand(fib_seed);

    shput_dup_value(config->config_data, "start_dx", "800.0");
    shput_dup_value(config->config_data, "start_dy", "800.0");
    shput_dup_value(config->config_data, "start_dz", "800.0");

    uint8_t size_code = 0;

    initialize_and_construct_grid(the_grid, POINT3D(204800, 204800, 204800));
    refine_grid(the_grid, 7);

    if(strcmp(scar_size, "big") == 0) {
        log_info("Loading Human Heart Edge with big scar\n");
        set_custom_mesh_with_bounds(the_grid, mesh_file, 2025252, 79100, 121000, 66700, 106000, 11200, 61400, "%lf,%lf,%lf,%lf,%d,%c\n");
        size_code = 0;
    } else if(strcmp(scar_size, "small") == 0) {
        log_info("Loading Human Heart Edge with small scar\n");
        set_custom_mesh_with_bounds(the_grid, mesh_file, 2025252, 30400, 81600, 59200, 103000, 13600, 48000, "%lf,%lf,%lf,%lf,%d,%c\n");
        size_code = 1;
    } else {
        log_error_and_exit("Function: initialize_grid_with_scar_edge, invalid scar size %s. Valid sizes are big or small. Exiting!\n", scar_size);
    }

    log_info("Cleaning grid\n");
    int i;
    for(i = 0; i < 7; i++) {
        derefine_grid_inactive_cells(the_grid);
    }

    refine_fibrotic_cells(the_grid);
    refine_fibrotic_cells(the_grid);
    refine_fibrotic_cells(the_grid);

    refine_border_zone_cells(the_grid);
    refine_border_zone_cells(the_grid);
    refine_border_zone_cells(the_grid);

    real_cpu scar_center_x;
    real_cpu scar_center_y;
    real_cpu scar_center_z;

    ////Fibrosis configuration

    // BIG SCAR
    if(size_code == 0) {
        scar_center_x = 95300.0;
        scar_center_y = 81600.0;
        scar_center_z = 36800.0;
    } else {
        scar_center_x = 52469.0;
        scar_center_y = 83225.0;
        scar_center_z = 24791.0;
    }

    real_cpu bz_size = 0.0;
    real_cpu dist;

    log_info("Using %u as seed\n", fib_seed);
    log_info("Calculating fibrosis using phi: %lf\n", phi);
    bool fibrotic, border_zone;

    FOR_EACH_CELL(the_grid) {
        if(cell->active) {
            fibrotic = FIBROTIC(cell);
            border_zone = BORDER_ZONE(cell);

            if(fibrotic) {
                cell->can_change = false;
                real_cpu p = (real_cpu)(rand()) / (RAND_MAX); // rand() has limited randomness
                if(p < phi)
                    cell->active = false;
            } else if(border_zone) {
                real_cpu center_x = cell->center.x;
                real_cpu center_y = cell->center.y;
                real_cpu center_z = cell->center.z;
                dist = sqrt((center_x - scar_center_x) * (center_x - scar_center_x) + (center_y - scar_center_y) * (center_y - scar_center_y) +
                            (center_z - scar_center_z) * (center_z - scar_center_z));
                if(dist > bz_size) {
                    bz_size = dist;
                }
            }
        }
    }

    FOR_EACH_CELL(the_grid) {
        if(cell->active) {
            border_zone = BORDER_ZONE(cell);
            if(border_zone) {
                real_cpu center_x = cell->center.x;
                real_cpu center_y = cell->center.y;
                real_cpu center_z = cell->center.z;
                dist = sqrt((center_x - scar_center_x) * (center_x - scar_center_x) + (center_y - scar_center_y) * (center_y - scar_center_y) +
                            (center_z - scar_center_z) * (center_z - scar_center_z));
                dist = dist / bz_size;

                real_cpu phi_local = phi - phi * dist;
                real_cpu p = (real_cpu)(rand()) / (RAND_MAX);

                if(p < phi_local) {
                    cell->active = false;
                }

                cell->can_change = false;
            }
        }
    }

    free(mesh_file);
    free(scar_size);

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_mouse_mesh) {

    char *mesh_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(mesh_file, config, "mesh_file");

    real_cpu start_h = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_h, config, "start_discretization");

    real_cpu maximum_discretization = start_h;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, maximum_discretization, config, "maximum_discretization");

    assert(the_grid);

    int num_volumes = 96195;

    log_info("Loading Mouse Heart Mesh\n");

    int num_loaded = set_custom_mesh_from_file(the_grid, mesh_file, num_volumes, start_h, 0, NULL);

    log_info("Read %d volumes from file: %s\n", num_loaded, mesh_file);

    free(mesh_file);

    if(start_h == 100.0) {
        //Do nothing
    }
    else if(start_h == 50.0) {
        log_info("Refining Mesh to 50um\n");
        refine_grid(the_grid, 1);
    } else if(start_h == 25.0) {
        log_info("Refining Mesh to 25um\n");
        refine_grid(the_grid, 2);
    } else if(start_h == 12.5) {
        log_info("Refining Mesh to 12.5um\n");
        refine_grid(the_grid, 3);
    } else {
        log_error("Invalid discretizations for this mesh. Valid discretizations are: 100um, 50um, 25um "
                  "or 12.5um. Using 100um!\n");
        start_h = 100.0;
    }

    the_grid->max_discretization = SAME_POINT3D(maximum_discretization);
    the_grid->start_discretization = SAME_POINT3D(start_h);

    return 1;
}

SET_CUSTOM_DATA_FOR_MESH(set_custom_data_for_atrial_percolation) {
    INITIALIZE_FIBROTIC_INFO(cell);
    FIBROTIC(cell) = custom_data[0] == 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_atrial_mesh_with_percolation) {

    char *mesh_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(mesh_file, config, "mesh_file");

    real_cpu start_h = 500.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, start_h, config, "original_discretization");

    real_cpu desired_h = 500.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, desired_h, config, "desired_discretization");

    int num_volumes = 1004741;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, num_volumes, config, "num_volumes");

    int num_loaded = set_custom_mesh_from_file(the_grid, mesh_file, num_volumes, start_h, 1, set_custom_data_for_atrial_percolation);

    log_info("Read %d volumes from file: %s\n", num_loaded, mesh_file);

    free(mesh_file);

    FOR_EACH_CELL(the_grid) {
        if(cell->active) {
            if(FIBROTIC(cell)) {
                cell->active = false;
            }
        }
    }

    the_grid->max_discretization = SAME_POINT3D(desired_h);
    the_grid->start_discretization = SAME_POINT3D(start_h);

    return num_loaded;
}

SET_CUSTOM_DATA_FOR_MESH(set_custom_data_for_scv_mesh) {

    INITIALIZE_SCV_INFO(cell);

    TISSUE_TYPE(cell) = custom_data[0];

    FIBROTIC(cell) = (custom_data[1] == 2);
    BORDER_ZONE(cell) = (custom_data[1] == 1);

    cell->sigma.fibers.f[0] = custom_data[2];
    cell->sigma.fibers.f[1] = custom_data[3];
    cell->sigma.fibers.f[2] = custom_data[4];
}

SET_CUSTOM_DATA_FOR_MESH(set_custom_data_for_hu_mesh) {
    INITIALIZE_FIBROTIC_INFO(cell);
    FIBROTIC(cell) = (custom_data[0] == 1);
}

//SET_CUSTOM_DATA_FOR_MESH(set_custom_data_for_hu_mesh_with_fibers) {
  //  INITIALIZE_FIBROTIC_INFO(cell);
   // FIBROTIC(cell) = (custom_data[0] == 1);

    //cell->sigma.fibers.f[0] = custom_data[1];
    //cell->sigma.fibers.f[1] = custom_data[2];
    //cell->sigma.fibers.f[2] = custom_data[3];
//}

SET_CUSTOM_DATA_FOR_MESH(set_custom_data_for_hu_mesh_with_fibers) {
    INITIALIZE_FIBROTIC_INFO(cell);
    // ordem antiga
    //FIBROTIC(cell) = (custom_data[0] == 1);

    //cell->sigma.fibers.f[0] = custom_data[1];
    //cell->sigma.fibers.f[1] = custom_data[2];
    //cell->sigma.fibers.f[2] = custom_data[3];
    
    FIBROTIC(cell) = (custom_data[4] == 1);

    cell->sigma.fibers.f[0] = custom_data[0];
    cell->sigma.fibers.f[1] = custom_data[1];
    cell->sigma.fibers.f[2] = custom_data[2];
    
    // informacao do u
    if (custom_data[3] < 0.45) {
      TISSUE_TYPE(cell) = 0; // ENDO
    }
    else if (custom_data[3] < 0.7) {
      TISSUE_TYPE(cell) = 1; // MID
    }
    else {
      TISSUE_TYPE(cell) = 2; // EPI
    }

    // TISSUE_TYPE(cell) = custom_data[3];
}

SET_SPATIAL_DOMAIN(initialize_grid_with_hu_mesh_with_scar) {

    char *mesh_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(mesh_file, config, "mesh_file");

    real_cpu start_h = 100.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, start_h, config, "original_discretization");

    real_cpu desired_h = 100.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, desired_h, config, "desired_discretization");

    assert(the_grid);

    log_info("Loading HU Heart Mesh with discretization: %lf\n", desired_h);

    uint32_t num_volumes = 514389;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(uint32_t, num_volumes, config, "num_volumes");

    uint32_t num_extra_fields = 1;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(uint32_t, num_extra_fields, config, "num_extra_fields");

    int num_loaded = 0;

    if(num_extra_fields == 1)
        num_loaded = set_custom_mesh_from_file(the_grid, mesh_file, num_volumes, start_h, num_extra_fields, set_custom_data_for_hu_mesh);
    else
        num_loaded = set_custom_mesh_from_file(the_grid, mesh_file, num_volumes, start_h, num_extra_fields, set_custom_data_for_hu_mesh_with_fibers);

    log_info("Read %d volumes from file (expected %d): %s\n", num_loaded, num_volumes, mesh_file);

    int num_refs_remaining = (int)(start_h / desired_h) - 1;

    if(num_refs_remaining > 0) {
        refine_grid(the_grid, num_refs_remaining);
    }

    the_grid->start_discretization = SAME_POINT3D(desired_h);

    free(mesh_file);

    return num_loaded;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_scv_mesh_with_scar_and_fibers) {

    char *mesh_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(mesh_file, config, "mesh_file");

    real_cpu start_h = 500.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, start_h, config, "original_discretization");

    real_cpu desired_h = 500.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, desired_h, config, "desired_discretization");

    assert(the_grid);

    log_info("Loading Atrial Heart Mesh with discretization: %lf\n", desired_h);

    uint32_t num_volumes = 514389;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(uint32_t, num_volumes, config, "num_volumes");

    uint32_t num_extra_fields = 5;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(uint32_t, num_extra_fields, config, "num_extra_fields");

    int num_loaded = set_custom_mesh_from_file(the_grid, mesh_file, num_volumes, start_h, num_extra_fields, set_custom_data_for_scv_mesh);

    log_info("Read %d volumes from file (expected %d): %s\n", num_loaded, num_volumes, mesh_file);

    int num_refs_remaining = (int)(start_h / desired_h) - 1;

    if(num_refs_remaining > 0) {
        refine_grid(the_grid, num_refs_remaining);
    }

    the_grid->start_discretization = SAME_POINT3D(desired_h);

    free(mesh_file);

    return num_loaded;
}

SET_SPATIAL_DOMAIN(set_perlin_square_mesh) {

    assert(the_grid);

    log_info("Loading perlin mesh\n");

    char *mesh_file = NULL;

    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(mesh_file, config, "mesh_file");

    real_cpu start_h = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, start_h, config, "start_discretization");

    real_cpu side_length = 0.0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, side_length, config, "side_length");

    size_t n_points = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(size_t, n_points, config, "n_points");

    set_cuboid_domain_mesh(the_grid, start_h, start_h, start_h, side_length, side_length, start_h);

    printf("Reading mesh file %s\n", mesh_file);
    set_custom_mesh(the_grid, mesh_file, n_points, "%lf,%lf,%lf,%lf,%d\n");

    FOR_EACH_CELL(the_grid) {
        if(cell->active) {
            if(FIBROTIC(cell)) {
                cell->active = false;
            }
        }
    }

    free(mesh_file);

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_with_plain_fibrotic_mesh_using_file) {

    set_square_mesh(config, the_grid);
    set_plain_fibrosis_using_file(the_grid, "fibrotic_positions.txt");

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_grid_hcm_mesh) {

    // TODO: we should put this in the grid data again
    //
    // the_grid -> start_discretization.x = SAME_POINT3D{start_discretization};
    //

    real_cpu discretization = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(real_cpu, discretization, config, "discretization");

    size_t size = 0;
    GET_PARAMETER_NUMERIC_VALUE_OR_REPORT_ERROR(size_t, size, config, "num_volumes");

    char *mesh_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(mesh_file, config, "mesh_file");

    // TODO: change this to the code above
    if(discretization == 500) {
        char *tmp = "500";
        shput_dup_value(config->config_data, "start_dx", tmp);
        shput_dup_value(config->config_data, "start_dy", tmp);
        shput_dup_value(config->config_data, "start_dz", tmp);
    } else if(discretization == 400) {
        char *tmp = "400";
        shput_dup_value(config->config_data, "start_dx", tmp);
        shput_dup_value(config->config_data, "start_dy", tmp);
        shput_dup_value(config->config_data, "start_dz", tmp);
    } else {
        log_error_and_exit("Discretization of %lf not supported for this mesh\n", discretization);
    }

    int n_steps;

    if(discretization == 500) {
        initialize_and_construct_grid(the_grid, POINT3D(128000, 128000, 128000));
        n_steps = 7;
    } else {
        initialize_and_construct_grid(the_grid, POINT3D(204800, 204800, 204800));
        n_steps = 8;
    }

    log_info("Refining the mesh\n");
    refine_grid(the_grid, n_steps);

    struct cell_node *grid_cell = the_grid->first_cell;
    FILE *file = fopen(mesh_file, "r");

    if(!file) {
        log_error_and_exit("Error opening mesh described in %s!!\n", mesh_file);
    }

    real_cpu **mesh_points = MALLOC_ARRAY_OF_TYPE(real_cpu *, size);

    for(int i = 0; i < size; i++) {
        mesh_points[i] = MALLOC_ARRAY_OF_TYPE(real_cpu, 4);
        if(mesh_points[i] == NULL) {
            log_error_and_exit("Failed to allocate memory\n");
        }
    }

    real_cpu maxy = 0.0;
    real_cpu maxz = 0.0;
    real_cpu miny = DBL_MAX;
    real_cpu minz = DBL_MAX;

    real_cpu *septum = MALLOC_ARRAY_OF_TYPE(real_cpu, size);
    real_cpu *tissue_type = MALLOC_ARRAY_OF_TYPE(real_cpu, size);

    uint32_t *element_id = MALLOC_ARRAY_OF_TYPE(uint32_t, size);
    uint32_t *node_id = MALLOC_ARRAY_OF_TYPE(uint32_t, size);

    real_cpu dummy;

    log_info("Setting mesh from file %s\n", mesh_file);
    int i = 0;
    while(i < size) {

        fscanf(file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d\n", &mesh_points[i][0], &mesh_points[i][1], &mesh_points[i][2], &dummy, &dummy, &dummy,
               &tissue_type[i], &septum[i], &element_id[i], &node_id[i]);

        // this is needed because the array mesh_points is sorted after reading the mesh file.
        mesh_points[i][3] = i;

        if(mesh_points[i][1] > maxy)
            maxy = mesh_points[i][1];
        if(mesh_points[i][2] > maxz)
            maxz = mesh_points[i][2];
        if(mesh_points[i][1] < miny)
            miny = mesh_points[i][1];
        if(mesh_points[i][2] < minz)
            minz = mesh_points[i][2];

        i++;
    }

    sort_vector(mesh_points, size); // we need to sort because inside_mesh perform a binary search

    real_cpu maxx = mesh_points[size - 1][0];
    real_cpu minx = mesh_points[0][0];
    int index;

    real_cpu x, y, z;
    while(grid_cell != 0) {
        x = grid_cell->center.x;
        y = grid_cell->center.y;
        z = grid_cell->center.z;

        if(x > maxx || y > maxy || z > maxz || x < minx || y < miny || z < minz) {
            grid_cell->active = false;
        } else {
            index = inside_mesh(mesh_points, x, y, z, 0, size - 1);

            if(index != -1) {
                grid_cell->active = true;
                int old_index = (int)mesh_points[index][3];

                INITIALIZE_HCM_INFO(grid_cell);

                HCM_SEPTUM(grid_cell) = (septum[old_index] == 2);
                HCM_TISSUE_TYPE(grid_cell) = tissue_type[old_index];
                HCM_ELEMENT_ID(grid_cell) = element_id[old_index];
                HCM_NODE_ID(grid_cell) = node_id[old_index];
                grid_cell->original_position_in_file = old_index;

            } else {
                grid_cell->active = false;
            }
        }
        grid_cell = grid_cell->next;
    }

    fclose(file);

    // deallocate memory
    for(int l = 0; l < size; l++) {
        free(mesh_points[l]);
    }

    free(mesh_points);
    free(mesh_file);

    the_grid->mesh_side_length.x = maxx;
    the_grid->mesh_side_length.y = maxy;
    the_grid->mesh_side_length.z = maxz;

    log_info("Cleaning grid\n");

    for(i = 0; i < n_steps; i++) {
        derefine_grid_inactive_cells(the_grid);
    }

    return 1;
}

SET_SPATIAL_DOMAIN(initialize_biatrial_mesh) {
        char *mesh_file = NULL;
        GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(mesh_file, config, "mesh_file");

        real_cpu start_h = 400.0;
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, start_h, config, "original_discretization");

        real_cpu desired_h = 400.0;
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(real_cpu, desired_h, config, "desired_discretization");

        assert(the_grid);

        log_info("Loading Atrial Heart Mesh with discretization: %lf\n", desired_h);

        uint32_t num_volumes = 1004741;
        GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(uint32_t, num_volumes, config, "num_volumes");

        int num_loaded = set_custom_mesh_from_file(the_grid, mesh_file, num_volumes, start_h, 0, NULL);

        log_info("Read %d volumes from file: %s\n", num_loaded, mesh_file);

        int num_refs_remaining = calc_num_refs(start_h, desired_h);

        if(num_refs_remaining > 0)
                refine_grid(the_grid, num_refs_remaining);

        the_grid->start_discretization = SAME_POINT3D(desired_h);

        free(mesh_file);

        return num_loaded;
}

SET_CUSTOM_DATA_FOR_MESH(set_custom_data_for_james_cube) {
    INITIALIZE_FIBROTIC_INFO(cell);
    FIBROTIC(cell) = (custom_data[4] == 1);
}

SET_SPATIAL_DOMAIN(initialize_grid_with_james_cube) {

    char *mesh_file = NULL;
    GET_PARAMETER_STRING_VALUE_OR_REPORT_ERROR(mesh_file, config, "mesh_file");

    int num_volumes = 1000000;
    GET_PARAMETER_NUMERIC_VALUE_OR_USE_DEFAULT(int, num_volumes, config, "num_volumes");

     real_cpu start_h = 300.0;
    int num_loaded = set_custom_mesh_from_file(the_grid, mesh_file, num_volumes, start_h, 5, set_custom_data_for_james_cube);

    log_info("Read %d volumes from file: %s\n", num_loaded, mesh_file);

    free(mesh_file);

    FOR_EACH_CELL(the_grid) {
        if(cell->active) {
            if(FIBROTIC(cell)) {
                cell->active = false;
            }
        }
    }

    the_grid->max_discretization = SAME_POINT3D(start_h);
    the_grid->start_discretization = SAME_POINT3D(start_h);

    return num_loaded;
}

SET_SPATIAL_DOMAIN(tt3_mixed) {
    return 1; //TODO: create a better example
}
