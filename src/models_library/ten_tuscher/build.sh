############## tt3_mixed ##############################
MODEL_FILE_CPU="tt3_mixed.c"
MODEL_FILE_GPU="tt3_mixed.cu"
COMMON_HEADERS="tt3_mixed.h"

COMPILE_MODEL_LIB "tt3_mixed" "$MODEL_FILE_CPU" "$MODEL_FILE_GPU" "$COMMON_HEADERS"
