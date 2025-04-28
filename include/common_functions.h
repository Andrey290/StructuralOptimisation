#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H

#include "structures.h"

// GENOTYPE FUNCTIONS
// genotype -> phenotype 
void genotype_phenotype(const double genotype[], double phenotype[]);
// phenotype -> genotype 
void phenotype_genotype(const double phenotype[], double genotype[]);
// initializers
void initialize_population(Population *pop);
void initialize_box(Box *box);
void initialize_room(Room *room);
void initialize_target(Target *target);
// functions needed for main circle in EntryPoint.c
void evaluate_population(Population *pop, Box *box, Room *room, Target *target, int k);
void select_parent(Population *pop, Individual *parent1, Individual *parent2);
void crossover(Individual *parent1, Individual *parent2, Individual *child);
void mutate(Individual *ind);
void replace_population(Population *pop, Population *new_pop);
void serialize_phenotype(char* buffer, const double phenotype[]);

#endif
