#include <unistd.h>

#include "../config/extra_data_config.h"
#include "../config_helpers/config_helpers.h"
#include "../libraries_common/common_data_structures.h"
#include "../domains_library/mesh_info_data.h"
#include "helper_functions.h"
#include "../logger/logger.h"
SET_EXTRA_DATA(set_extra_data_for_scv_mesh) {
    log_info("-----------------------------------xxxxxxx\n");

    uint32_t num_active_cells = the_grid->num_active_cells;

    *extra_data_size = sizeof(uint32_t)*(num_active_cells);

    uint32_t *mapping = (uint32_t*)malloc(*extra_data_size);

    struct cell_node ** ac = the_grid->active_cells;
    log_info("num_active_cells na SET_EXTRA_DATA = %d \n" ,num_active_cells);

    int i;

    OMP(parallel for)
    for (i = 0; i < num_active_cells; i++) {
        

        mapping[i] = TISSUE_TYPE(ac[i]); //endo 0; mid  1; epi  2
    }

    return (void*)mapping;
}

// Initial condition - 'libToRORd_fkatp_mixed_endo_mid_epi.so' + transmurality (cable and cuboid)
SET_EXTRA_DATA(set_extra_data_mixed_minimal_model_square) {

    uint32_t num_active_cells = the_grid->num_active_cells;
    real side_length = the_grid->mesh_side_length.x;
    struct cell_node ** ac = the_grid->active_cells;

    *extra_data_size = sizeof(uint32_t)*(num_active_cells);
    uint32_t *mapping = (uint32_t*)malloc(*extra_data_size);

    // The percentages were taken from the ToRORd paper (Transmural experiment)
    real side_length_endo = side_length*0.45;
    real side_length_mid = side_length_endo + side_length*0.25;
    real side_length_epi = side_length_mid + side_length*0.3;
    log_info("num_active_cells na SET_EXTRA_DATA = %d \n" ,num_active_cells);
    int i;

    OMP(parallel for)
    for (i = 0; i < num_active_cells; i++) {
        real center_x = ac[i]->center.x;
        // ENDO
        if (center_x < side_length_endo)
            mapping[i] = 0.0;
        // MID
        else if (center_x >= side_length_endo && center_x < side_length_mid)
            mapping[i] = 1.0;
        // EPI
        else
            mapping[i] = 2.0;
    }

    return (void*)mapping;
}