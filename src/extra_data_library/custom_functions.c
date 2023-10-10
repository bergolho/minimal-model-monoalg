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