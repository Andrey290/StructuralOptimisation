//---SELF_HEADER---//

#include "common_functions.h"

//---HEADERS---//

#include "simulation.h"
#include "config.h"
#include "structures.h"
#include "logger.h"

//---LIBRARIES---//
#include <stdlib.h>
#include <math.h>

// genotype -> phenotype
void genotype_phenotype(const double genotype[], double phenotype[]) {
    for (int i = 0; i < NUM_GENES; i++) {
        phenotype[i] = 0.0;
       	for (int j = 0; j < NUM_GENES; j++) {
	    phenotype[i] += genotype[j] * coefficients[i][j];
	}
    }
}

// phenotype -> genotype
void phenotype_genotype(const double phenotype[], double genotype[]) {
    for (int j = 0; j < NUM_GENES; j++) {
	double sum = 0.0;
        for (int i = 0; i < NUM_GENES; i++) {
           sum += inverse_coefficients[j][i] * phenotype[i];
	}
	genotype[j] = sum;
    }
}

// Population initialization
void initialize_population(Population *pop) {
    for (int i = 0; i < POPULATION_SIZE; i++) { 
        //-------------------------FITNESS------------------//
	pop->individuals[i].fitness = 0;
	//------------------------PHENOTYPE-----------------//
	pop->individuals[i].phenotype[0] = LINK_1_MIN + (double)rand() / RAND_MAX * (LINK_1_MAX - LINK_1_MIN); // Добавляем некую рандомную долю допустимого угла к 
	pop->individuals[i].phenotype[1] = LINK_1_MAX + (double)rand() / RAND_MAX * (LINK_1_MAX - LINK_1_MIN); // стартовым значениям ограничений
	pop->individuals[i].phenotype[2] = LINK_2_MIN + (double)rand() / RAND_MAX * (LINK_2_MAX - LINK_2_MIN);
	pop->individuals[i].phenotype[3] = LINK_2_MAX + (double)rand() / RAND_MAX * (LINK_2_MAX - LINK_2_MIN);
	pop->individuals[i].phenotype[4] = LINK_3_MIN + (double)rand() / RAND_MAX * (LINK_3_MAX - LINK_3_MIN);
	pop->individuals[i].phenotype[5] = LINK_3_MAX + (double)rand() / RAND_MAX * (LINK_3_MAX - LINK_3_MIN);

	pop->individuals[i].phenotype[6] = LINK_2_START_LEN + ((double)rand() / RAND_MAX) * LINK_2_START_LEN; // Добавляем к стандартной стартовой длине неокоторую долю 
	pop->individuals[i].phenotype[7] = LINK_3_START_LEN + ((double)rand() / RAND_MAX) * LINK_3_START_LEN; // стандартной стартовой длины.
	//-------------------------GENOTYPE-----------------//
        phenotype_genotype(pop->individuals[i].phenotype, pop->individuals[i].genotype); // here we calculate genotype (we do not need this foo any more)

        char phenotype_str[256];
        serialize_phenotype(phenotype_str, pop->individuals[i].phenotype);
        LOG_EVENT("INIT_PHENOTYPE", "GENERATION=%d,INDIVIDUAL=%d,PHENOTYPE=[%s]", pop->generation, i, phenotype_str);
    }
    // LOG //
    LOG_EVENT("INIT_POPULATION", "GENERATION=%d", pop->generation);
}

// Box initialization
void initialize_box(Box *box) {
    // Botoom
    // // plane
    // // // normal
    box->walls[0].plane.normal[0] =  0;
    box->walls[0].plane.normal[1] =  0;
    box->walls[0].plane.normal[2] = -1;
    // // // distance
    box->walls[0].plane.distance = BOX_Z_COORD + -BOX_DIMENSION / 2;
    // // min bonds
    box->walls[0].min_bounds[0] = BOX_X_COORD + -BOX_DIMENSION / 2;  // x
    box->walls[0].min_bounds[1] = BOX_Y_COORD + -BOX_DIMENSION / 2;  // y
    box->walls[0].min_bounds[2] = BOX_Z_COORD + -BOX_DIMENSION / 2;  // z
    // // max bonds
    box->walls[0].max_bounds[0] = BOX_X_COORD +  BOX_DIMENSION / 2;  // x
    box->walls[0].max_bounds[1] = BOX_Y_COORD +  BOX_DIMENSION / 2;  // y
    box->walls[0].max_bounds[2] = BOX_Z_COORD + -BOX_DIMENSION / 2;  // z
    
    // West wall
    // // plane
    // // // normal
    box->walls[1].plane.normal[0] =  0;
    box->walls[1].plane.normal[1] = -1;
    box->walls[1].plane.normal[2] =  0;
    // // // distance
    box->walls[1].plane.distance = BOX_Y_COORD + -BOX_DIMENSION / 2;
    // // min bonds
    box->walls[1].min_bounds[0] = BOX_X_COORD + -BOX_DIMENSION / 2;  // x
    box->walls[1].min_bounds[1] = BOX_Y_COORD + -BOX_DIMENSION / 2;  // y
    box->walls[1].min_bounds[2] = BOX_Z_COORD + -BOX_DIMENSION / 2;  // z
    // // max bonds
    box->walls[1].max_bounds[0] = BOX_X_COORD +  BOX_DIMENSION / 2;  // x
    box->walls[1].max_bounds[1] = BOX_Y_COORD + -BOX_DIMENSION / 2;  // y
    box->walls[1].max_bounds[2] = BOX_Z_COORD +  BOX_DIMENSION / 2;  // z
    
    // South wall
    // // plane
    // // // normal
    box->walls[2].plane.normal[0] =  1;
    box->walls[2].plane.normal[1] =  0;
    box->walls[2].plane.normal[2] =  0;
    // // // distance
    box->walls[2].plane.distance = BOX_X_COORD +  BOX_DIMENSION / 2;
    // // min bonds
    box->walls[2].min_bounds[0] = BOX_X_COORD +  BOX_DIMENSION / 2;  // x
    box->walls[2].min_bounds[1] = BOX_Y_COORD + -BOX_DIMENSION / 2;  // y
    box->walls[2].min_bounds[2] = BOX_Z_COORD + -BOX_DIMENSION / 2;  // z
    // // max bonds
    box->walls[2].max_bounds[0] = BOX_X_COORD +  BOX_DIMENSION / 2;  // x
    box->walls[2].max_bounds[1] = BOX_Y_COORD +  BOX_DIMENSION / 2;  // y
    box->walls[2].max_bounds[2] = BOX_Z_COORD +  BOX_DIMENSION / 2;  // z

    // East wall
    // // plane
    // // // normal
    box->walls[3].plane.normal[0] =  0;
    box->walls[3].plane.normal[1] =  1;
    box->walls[3].plane.normal[2] =  0;
    // // // distance
    box->walls[3].plane.distance = BOX_Y_COORD +  BOX_DIMENSION / 2;
    // // min bonds
    box->walls[3].min_bounds[0] = BOX_X_COORD + -BOX_DIMENSION / 2;  // x
    box->walls[3].min_bounds[1] = BOX_Y_COORD +  BOX_DIMENSION / 2;  // y
    box->walls[3].min_bounds[2] = BOX_Z_COORD + -BOX_DIMENSION / 2;  // z
    // // max bonds
    box->walls[3].max_bounds[0] = BOX_X_COORD +  BOX_DIMENSION / 2;  // x
    box->walls[3].max_bounds[1] = BOX_Y_COORD +  BOX_DIMENSION / 2;  // y
    box->walls[3].max_bounds[2] = BOX_Z_COORD +  BOX_DIMENSION / 2;  // z

    // North wall
    // // plane
    // // // normal
    box->walls[4].plane.normal[0] = -1;
    box->walls[4].plane.normal[1] =  0;
    box->walls[4].plane.normal[2] =  0;
    // // // distance
    box->walls[4].plane.distance = BOX_X_COORD + -BOX_DIMENSION / 2;
    // // min bonds
    box->walls[4].min_bounds[0] = BOX_X_COORD + -BOX_DIMENSION / 2;  // x
    box->walls[4].min_bounds[1] = BOX_Y_COORD + -BOX_DIMENSION / 2;  // y
    box->walls[4].min_bounds[2] = BOX_Z_COORD + -BOX_DIMENSION / 2;  // z
    // // max bonds
    box->walls[4].max_bounds[0] = BOX_X_COORD + -BOX_DIMENSION / 2;  // x
    box->walls[4].max_bounds[1] = BOX_Y_COORD +  BOX_DIMENSION / 2;  // y
    box->walls[4].max_bounds[2] = BOX_Z_COORD +  BOX_DIMENSION / 2;  // z
}

// Room initialization
void initialize_room(Room *room) {
    // Botoom
    // // plane
    // // // normal
    room->walls[0].plane.normal[0] =  0;
    room->walls[0].plane.normal[1] =  0;
    room->walls[0].plane.normal[2] = -1;
    // // // distance
    room->walls[0].plane.distance = ROOM_Z_COORD + -ROOM_HEIGHT;
    // // min bonds
    room->walls[0].min_bounds[0] = ROOM_X_COORD + -ROOM_DIMENSION / 2;  // x
    room->walls[0].min_bounds[1] = ROOM_Y_COORD + -ROOM_DIMENSION / 2;  // y
    room->walls[0].min_bounds[2] = ROOM_Z_COORD + -ROOM_HEIGHT    / 2;  // z
    // // max bonds
    room->walls[0].max_bounds[0] = ROOM_X_COORD +  ROOM_DIMENSION / 2;  // x
    room->walls[0].max_bounds[1] = ROOM_Y_COORD +  ROOM_DIMENSION / 2;  // y
    room->walls[0].max_bounds[2] = ROOM_Z_COORD + -ROOM_HEIGHT    / 2;  // z
    
    // West wall
    // // plane
    // // // normal
    room->walls[1].plane.normal[0] =  0;
    room->walls[1].plane.normal[1] = -1;
    room->walls[1].plane.normal[2] =  0;
    // // // distance
    room->walls[1].plane.distance = ROOM_Y_COORD + -ROOM_DIMENSION / 2;
    // // min bonds
    room->walls[1].min_bounds[0] = ROOM_X_COORD + -ROOM_DIMENSION / 2;  // x
    room->walls[1].min_bounds[1] = ROOM_Y_COORD + -ROOM_DIMENSION / 2;  // y
    room->walls[1].min_bounds[2] = ROOM_Z_COORD + -ROOM_HEIGHT    / 2;  // z
    // // max bonds
    room->walls[1].max_bounds[0] = ROOM_X_COORD +  ROOM_DIMENSION / 2;  // x
    room->walls[1].max_bounds[1] = ROOM_Y_COORD + -ROOM_DIMENSION / 2;  // y
    room->walls[1].max_bounds[2] = ROOM_Z_COORD +  ROOM_HEIGHT    / 2;  // z
    
    // South wall
    // // plane
    // // // normal
    room->walls[2].plane.normal[0] =  1;
    room->walls[2].plane.normal[1] =  0;
    room->walls[2].plane.normal[2] =  0;
    // // // distance
    room->walls[2].plane.distance = ROOM_X_COORD +  ROOM_DIMENSION / 2;
    // // min bonds
    room->walls[2].min_bounds[0] = ROOM_X_COORD +  ROOM_DIMENSION / 2;  // x
    room->walls[2].min_bounds[1] = ROOM_Y_COORD + -ROOM_DIMENSION / 2;  // y
    room->walls[2].min_bounds[2] = ROOM_Z_COORD + -ROOM_HEIGHT    / 2;  // z
    // // max bonds
    room->walls[2].max_bounds[0] = ROOM_X_COORD +  ROOM_DIMENSION / 2;  // x
    room->walls[2].max_bounds[1] = ROOM_Y_COORD +  ROOM_DIMENSION / 2;  // y
    room->walls[2].max_bounds[2] = ROOM_Z_COORD +  ROOM_HEIGHT    / 2;  // z

    // East wall
    // // plane
    // // // normal
    room->walls[3].plane.normal[0] =  0;
    room->walls[3].plane.normal[1] =  1;
    room->walls[3].plane.normal[2] =  0;
    // // // distance
    room->walls[3].plane.distance = ROOM_Y_COORD +  ROOM_DIMENSION / 2;
    // // min bonds
    room->walls[3].min_bounds[0] = ROOM_X_COORD + -ROOM_DIMENSION / 2;  // x
    room->walls[3].min_bounds[1] = ROOM_Y_COORD +  ROOM_DIMENSION / 2;  // y
    room->walls[3].min_bounds[2] = ROOM_Z_COORD + -ROOM_HEIGHT    / 2;  // z
    // // max bonds
    room->walls[3].max_bounds[0] = ROOM_X_COORD +  ROOM_DIMENSION / 2;  // x
    room->walls[3].max_bounds[1] = ROOM_Y_COORD +  ROOM_DIMENSION / 2;  // y
    room->walls[3].max_bounds[2] = ROOM_Z_COORD +  ROOM_HEIGHT    / 2;  // z

    // North wall
    // // plane
    // // // normal
    room->walls[4].plane.normal[0] = -1;
    room->walls[4].plane.normal[1] =  0;
    room->walls[4].plane.normal[2] =  0;
    // // // distance
    room->walls[4].plane.distance = ROOM_X_COORD + -ROOM_DIMENSION / 2;
    // // min bonds
    room->walls[4].min_bounds[0] = ROOM_X_COORD + -ROOM_DIMENSION / 2;  // x
    room->walls[4].min_bounds[1] = ROOM_Y_COORD + -ROOM_DIMENSION / 2;  // y
    room->walls[4].min_bounds[2] = ROOM_Z_COORD + -ROOM_HEIGHT    / 2;  // z
    // // max bonds
    room->walls[4].max_bounds[0] = ROOM_X_COORD + -ROOM_DIMENSION / 2;  // x
    room->walls[4].max_bounds[1] = ROOM_Y_COORD +  ROOM_DIMENSION / 2;  // y
    room->walls[4].max_bounds[2] = ROOM_Z_COORD +  ROOM_HEIGHT    / 2;  // z

    // Cap
    // // plane
    // // // normal
    room->walls[5].plane.normal[0] =  0;
    room->walls[5].plane.normal[1] =  0;
    room->walls[5].plane.normal[2] =  1;
    // // // distance
    room->walls[5].plane.distance = ROOM_Z_COORD +  ROOM_HEIGHT / 2;
    // // min bonds
    room->walls[5].min_bounds[0] = ROOM_X_COORD + -ROOM_DIMENSION / 2;  // x
    room->walls[5].min_bounds[1] = ROOM_Y_COORD + -ROOM_DIMENSION / 2;  // y
    room->walls[5].min_bounds[2] = ROOM_Z_COORD +  ROOM_HEIGHT    / 2;  // z
    // // max bonds
    room->walls[5].max_bounds[0] = ROOM_X_COORD +  ROOM_DIMENSION / 2;  // x
    room->walls[5].max_bounds[1] = ROOM_Y_COORD +  ROOM_DIMENSION / 2;  // y
    room->walls[5].max_bounds[2] = ROOM_Z_COORD +  ROOM_HEIGHT    / 2;  // z
}

// Target initialization
void initialize_target(Target *target) {
    target->position[0] = TARGET_X;
    target->position[1] = TARGET_Y;
    target->position[2] = TARGET_Z;

    target->axis[0] = 0;
    target->axis[1] = 0;
    target->axis[2] = 1;
}

void evaluate_population(Population *pop, Box *box, Room *room, Target *target, int k) {
    // for every pop->individuals 
    for (int i = 0; i < POPULATION_SIZE; i++) {
        // angle borders calculation
	// ADD FURTHER
       
        // length_metric calculation
        double length_metric = 1 / (pop->individuals[i].phenotype[6] + pop->individuals[i].phenotype[7]);
    
        // SIMULATION //
	// precision_metric calculation
        double precision_metric = 0.0;
        int collision_count = 0;
	simulation(&precision_metric, &collision_count, &(pop->individuals[i]), box, room, target, i); // distance and service angle
	double penalty = collision_count * 1000;

	double antiprecision_metric = 1000 / precision_metric;
        // OUTPUT //
        pop->individuals[i].fitness = length_metric + antiprecision_metric - penalty;
	// LOG //
	char phenotype_str[256];
        serialize_phenotype(phenotype_str, pop->individuals[i].phenotype);
        LOG_EVENT("FITNESS", "INDIVIDUAL=%d,FITNESS=%.3f,PRECISION=%.3f,COLLISIONS=%d,PHENOTYPE=[%s]", 
			i, pop->individuals[i].fitness, precision_metric, collision_count, phenotype_str);  	
    }
}

void select_parent(Population *pop, Individual *parent1, Individual *parent2) {
    // TOURNAMENT SELECTION
    // турнирная сселекция может упустить лучшего!
    int tournament_size =  POPULATION_SIZE; // P_S / 2
    
    // first ancestor (parent1) selection
    int first_best_idx = rand() % POPULATION_SIZE;
    for (int i = 0; i < tournament_size; i++) {
        
	int contender_idx = rand() % POPULATION_SIZE;
	if (pop->individuals[contender_idx].fitness > pop->individuals[first_best_idx].fitness) {
            first_best_idx = contender_idx;
	}
	
	/*
	if (pop->individuals[i].fitness > -1000 && pop->individuals[i].fitness > pop->individuals[first_best_idx].fitness) {
            first_best_idx = i;
	}
	*/

    }
    *parent1 = pop->individuals[first_best_idx];

    // second ancestor (parent2) selection
    int second_best_idx = rand() % POPULATION_SIZE;
    for (int i = 1; i < tournament_size; i++) {
            
	    int contender_idx = rand() % POPULATION_SIZE;
	    if (pop->individuals[contender_idx].fitness > pop->individuals[second_best_idx].fitness && contender_idx != first_best_idx) {
                second_best_idx = contender_idx;
	    }
	    
	    /*
            if (pop->individuals[i].fitness > -1000 && pop->individuals[i].fitness > pop->individuals[second_best_idx].fitness && i != first_best_idx) {
                second_best_idx = i;
	    }
	    */
    }
    *parent2 = pop->individuals[second_best_idx];
    // LOG //
    LOG_EVENT("SELECT_PARENT",
    "PARENT1=%d,PARENT2=%d", first_best_idx, second_best_idx
    );
}

    /* // ELITISM SELECTION (ADD THIS IN ALT FORK)
    Individual best = pop->individuals[0];
    Individual prev_best;
    Individual contender;
    for (int i = 1; i < POPULATION_SIZE; i++) {
        contender = pop->individuals[i];
        if (contender.fitness > best.fitness) {
            prev_best = best;
	    best = contender;
        }
    }
    parent1 = best;
    parent2 = prev_best;*/


void crossover(Individual *parent1, Individual *parent2, Individual *child) {
    int HALF_GENES = NUM_GENES / 2;
    for (int i = 0; i < HALF_GENES; i++) {
	child->genotype[i] = parent1->genotype[i];
    }
    for (int i = HALF_GENES; i < NUM_GENES; i++) {
	child->genotype[i] = parent2->genotype[i];
    }
    // DECODING GENOTYPE RIGHT HERE
    genotype_phenotype(child->genotype, child->phenotype);
}

void mutate(Individual *ind) {
    // MUTATE SINGLE RANDOM GEN
    //float sign = (rand() % 2) ? 1.0f : -1.0f;
    int rand_gen = rand() % NUM_GENES; // choose single gen
    //ind->genotype[rand_gen] += sign * 0.1 * ind->genotype[rand_gen];
    float mut_rate = 0.1f;
    float mutation = mut_rate * ind->genotype[rand_gen];
    ind->genotype[rand_gen] += (rand() % 2) ? mutation : -mutation; // рандомно мутируем или в плюс или в миус
}

void replace_population(Population *pop, Population *new_pop) {
    int prev_generation = pop->generation;
    *pop = *new_pop;
    pop->generation = prev_generation + 1;
    // LOG //
    LOG_EVENT("REPLACE_POPULATION", "GENERATION=%d", pop->generation);
}

// Сериализация фенотипа в строку
void serialize_phenotype(char* buffer, const double phenotype[]) {
    int offset = 0;
    for (int i = 0; i < NUM_GENES; i++) {
        offset += snprintf(buffer + offset, 256 - offset, "%.3f,", phenotype[i]);
    }
    buffer[offset - 1] = '\0'; // Удаляем последнюю запятую
}
