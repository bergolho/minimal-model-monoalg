############## torord_extra ##############################
MODEL_FILE_CPU="torord_extra.c"
MODEL_FILE_GPU="torord_extra.cu"
COMMON_HEADERS="torord_extra.h"

COMPILE_MODEL_LIB "torord_extra" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"
