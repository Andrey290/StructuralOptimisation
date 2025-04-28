#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "config.h"
#include <stdbool.h>

//---MANIPULATORS---//

typedef struct {
	double fitness;              // Оценка, определяемая функцией оценивания
	double genotype[NUM_GENES];  // В файле конфигов надо определён матрицы генотип->фенотип и фенотип->генотип
	double phenotype[NUM_GENES]; // Мы будем использовать генотип той же длины, что и фенотип, но при этом каждый ген будет иметь влияние на несколько фенов
} Individual;

typedef struct {
	Individual individuals[POPULATION_SIZE];
	int generation;
} Population;

typedef struct {
    double start[3];
    double end[3];
} Segment;

typedef struct {
    double generalized_coordinates[LINKS];
    double quality;
} Pose;

//---ENVIROMENT---//

// Plane
typedef struct {
    double normal[3];
    double distance; // from (0, 0, 0)
} Plane;

// Wall
typedef struct {
    Plane plane;
    double min_bounds[3]; // левый нижний угол
    double max_bounds[3]; // правый верхний угол
} Wall;

// Box
typedef struct {
    Wall walls[5]; // 5 box walls
} Box;

// Room
typedef struct {
    Wall walls[6]; // 6 room walls
} Room;

// Target
typedef struct {
    double position[3];
    double axis[3];
} Target;

#endif
