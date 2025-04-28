#ifndef SIMULATION_H
#define SIMULATION_H

#include "structures.h"

void simulation(double *precision_metric, int *collision_count, const Individual *ind, const Box *box, const Room *room, const Target *target, int id);

int is_point_in_wall(const double point[3], const Wall *wall);

int check_wall_collision(const double point[3], const Wall *wall);

void copy_matrixes(double v[4][4], const double c[4][4]);
	
void mult_matrixes(double result[4][4], const double a[4][4], const double b[4][4]);

void initialize_pose(Pose *pose, const Individual *ind);

void copy_pose(Pose *new_pose, const Pose *pose);

int dh_parametrs_filling(double **array, const Individual *ind, double *generalized_coordinates);

int collision_counter(const double beg[3], const double end[3], const Box *box, const Room *room);

void init_dh_matrix(double A[4][4], double **array, int j);


#endif


