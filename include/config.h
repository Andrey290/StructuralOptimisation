#ifndef CONFIG_H
#define CONFIG_H

//////////////////////
// GENETIC FEATURES //
//////////////////////

#define NUM_GENES 8
#define BASE 100

extern const double          coefficients[NUM_GENES][NUM_GENES];
extern const double inverse_coefficients [NUM_GENES][NUM_GENES];

#define POPULATION_SIZE 300
#define MAX_GENERATIONS 50
#define MUTATION_RATE 0.01
#define CROSSOVER_RATE 0.7

/////////////////////////
// MECHANICAL FEATURES //
/////////////////////////

#define LINKS 7                
#define DH_PARAMS_COLS 4 // !

//*** Links ***
#define LINK_0_LEN 1          // Fixed
#define LINK_1_LEN 1          // Fixed
#define LINK_2_START_LEN 5    // Init
#define LINK_3_START_LEN 5    // Init
#define GRAB_LEN   1          // Fixed

//*** Limits ***
#define LINK_0_MIN    -180     // Fixed
#define LINK_0_MAX     180     // Fixed
		        
#define LINK_1_MIN    -60      // Init
#define LINK_1_MAX     60      // Init

#define LINK_2_MIN     0       // Init
#define LINK_2_MAX     90      // Init

#define LINK_3_MIN    -90      // Init
#define LINK_3_MAX     0       // Init

#define GRAB_NUT_MIN  -30      // Fixed
#define GRAB_NUT_MAX   30      // Fixed
#define GRAB_PRC_MIN  -30      // Fixed
#define GRAB_PRC_MAX   30      // Fixed
#define GRAB_SLF_MIN  -180      // Fixed
#define GRAB_SLF_MAX   180      // Fixed

//////////////
// ROOM&BOX //
//////////////

#define BOX_DIMENSION 6
#define BOX_X_COORD 8
#define BOX_Y_COORD 13
#define BOX_Z_COORD 3

#define ROOM_DIMENSION 16
#define ROOM_HEIGHT 8
#define ROOM_X_COORD 8
#define ROOM_Y_COORD 8
#define ROOM_Z_COORD 4

#define TARGET_X 8
#define TARGET_Y 13
#define TARGET_Z 1

/////////////////////
// NEWTON FEATURES //
/////////////////////

#define MAX_ITERATIONS 10
#define TOLERANCE 1e-6

/////////////////////////
// SIMULATION FEATURES //
/////////////////////////

#define STEP_DEGREES  5 // Degrees of step

//////////////
// TCHNICAL //
//////////////

static const double ZERO_ARRAY[LINKS] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

#endif

