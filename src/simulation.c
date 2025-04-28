//---SELF_HEADER---//

#include "simulation.h"

//---HEADERS---//

#include "config.h"
#include "structures.h"
#include "logger.h"

//---LIBRARIES---//
#include <stdlib.h>
#include <math.h>
 
void simulation(double *precision_metric, int *collision_count, const Individual *ind, const Box *box, const Room *room, const Target *target, int id) {
    // LOG //
    LOG_EVENT("SIMULATION_START", "INDIVIDUAL=%d", id);
    //-------ALLOCATION-------//
    // DH parametrs
    int rows = LINKS;
    int cols = DH_PARAMS_COLS; 
    double **array = malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++) {
	array[i] = malloc(cols * sizeof(double));
    }
    
    //-------CALCULATIONS-------//
    
    // Градиентный спуск для оптимизации позы
    const int max_iterations = 100;
    double learning_rate = 1;
    double epsilon = 1e-3; // Приращение
    
    Pose pose;
    initialize_pose(&pose, ind); 
    
    // |||| Этот цикл приводит позу pose в наилучшее по евклидовой метрике положение относительно цели |||| //
    // Важно добавить в будущем оценку не только евклидова расстояния, но и угла между осями схвата и цели. !!!
    for (int iter = 0; iter < max_iterations; iter++) {
        double gradient[LINKS] = {0};

        // || Этот цикл вычисляет градиент на данном шаге приближения || //
        for (int i = 0; i < LINKS; i++) {
            // Сохраняем изменённое значение
            Pose temp_pose;
            copy_pose(&temp_pose, &pose);
	    temp_pose.generalized_coordinates[i] += epsilon;
            dh_parametrs_filling(array, ind, temp_pose.generalized_coordinates);  // Обновление матриц DH
        
            // Пересчет качества с новыми параметрами
            double T[4][4] = {
		    {1,0,0,8},
		    {0,1,0,3},
		    {0,0,1,0},
		    {0,0,0,1}
	    };
            
	    for (int j = 0; j < LINKS; j++) {
                double A[4][4];               // Init single transformation matrix for this link
                init_dh_matrix(A, array, j);  // Fill it by values for this link from changed DH matrix

                double Tnew[4][4];            
                mult_matrixes(Tnew, (const double(*)[4])T, (const double(*)[4])A);
                copy_matrixes((double(*)[4])T, (const double(*)[4])Tnew);
            }
            
	    // We received effector's coords --> now computing quality increment
            double dx = T[0][3] - target->position[0];
            double dy = T[1][3] - target->position[1];
            double dz = T[2][3] - target->position[2];

            double new_quality = sqrt(dx*dx + dy*dy + dz*dz);
        
            // Computing derivative
            gradient[i] = (new_quality - pose.quality) / epsilon;
	    // если новое расстояние до цели уменьшилось по сравнению со текущим, то производая отрицательная
	    // так как производная показывает направление роста значения функции, то нам надо идти против производной
        }
	// Итак, мы вычислили градиент на этом шаге приближения
        // Теперь нам надо обновить обобщённые угловые координаты нашей рабочей позы
	// || Этот цикл применяет производные умноженные на шаг к текущим угловым координатам || //
	// Вообще гооворя, шаг learning_rate надо менять !!!
        for (int i = 0; i < LINKS; i++) {
            pose.generalized_coordinates[i] -= learning_rate * gradient[i];
        
            // Применение ограничений (для трёх координат ограничения берутся не из конфигов, а из индивида)
            switch(i) {
                case 0:
                    pose.generalized_coordinates[i] = fmax(LINK_0_MIN, fmin(pose.generalized_coordinates[i], LINK_0_MAX));
                    break;
                case 1:
                    pose.generalized_coordinates[i] = fmax(ind->phenotype[0], fmin(pose.generalized_coordinates[i], ind->phenotype[1]));
                    break;
                case 2:
                    pose.generalized_coordinates[i] = fmax(ind->phenotype[2], fmin(pose.generalized_coordinates[i], ind->phenotype[3]));
                    break;
		case 3:
                    pose.generalized_coordinates[i] = fmax(ind->phenotype[4], fmin(pose.generalized_coordinates[i], ind->phenotype[5]));
                    break;
                case 4:
                    pose.generalized_coordinates[i] = fmax(GRAB_NUT_MIN, fmin(pose.generalized_coordinates[i], GRAB_NUT_MAX));
                    break;
                case 5:
                    pose.generalized_coordinates[i] = fmax(GRAB_PRC_MIN, fmin(pose.generalized_coordinates[i], GRAB_PRC_MAX));
                    break;
                case 6:
                    pose.generalized_coordinates[i] = fmax(GRAB_SLF_MIN, fmin(pose.generalized_coordinates[i], GRAB_SLF_MAX));
                    break;
            }
        }

        // После произведённых изменений надо ещё раз пересчитать quality позы pose
	// И если это 10-я итерация, то посчитать число коллизий
        ///////////////////////////////////////////////////////////////////////////
	dh_parametrs_filling(array, ind, pose.generalized_coordinates);  // Обновление матриц DH
        double T[4][4] = {
            {1,0,0,8},
	    {0,1,0,3},
	    {0,0,1,0},
	    {0,0,0,1}
        };
	for (int j = 0; j < LINKS; j++) {
            double A[4][4];               // Init single transformation matrix for this link
            init_dh_matrix(A, array, j);  // Fill it by values for this link from changed DH matrix
            double Tnew[4][4];            
            mult_matrixes(Tnew, (const double(*)[4])T, (const double(*)[4])A);
	    // Если это финальная итерация
	    if (iter == max_iterations - 1) {
		double beg[3] = {   T[0][3],    T[1][3],    T[2][3]};
		double end[3] = {Tnew[0][3], Tnew[1][3], Tnew[2][3]};
		// проверяем, сколько плоскостей пересекает одно звено
                (*collision_count) += collision_counter(beg, end, box, room);
		// LOG //
		LOG_EVENT("COLLISION_CHECK",
                "SEGMENT=(%.2f,%.2f,%.2f)->(%.2f,%.2f,%.2f),COUNT=%d",
                beg[0], beg[1], beg[2], end[0], end[1], end[2], collision_count
                );
	    }
            copy_matrixes((double(*)[4])T, (const double(*)[4])Tnew);
        }
        double dx = T[0][3] - target->position[0];
        double dy = T[1][3] - target->position[1];
        double dz = T[2][3] - target->position[2];
	pose.quality = sqrt(dx*dx + dy*dy + dz*dz);
        *precision_metric = pose.quality;
        ///////////////////////////////////////////////////////////////////////////

	// А здесь мы как раз меняем шаг, но коэффициент уменьшения "магический". Его надо обосновать
        learning_rate /= iter;

	// LOG //
	LOG_EVENT("GRADIENT_STEP",
        "ITER=%d,QUALITY=%.3f,LEARNING_RATE=%.3f",
        iter, pose.quality, learning_rate
        );
    }

    //-------DEALLOCATION-------//
    for (int i = 0; i < rows; i++) {
	free(array[i]);
    }
    free(array);
}

//---COLLISIONS---//

// Проверка, лежит ли точка внутри прямоугольной стены (в её локальной системе координат)
int is_point_in_wall(const double point[3], const Wall *wall) {
    for (int i = 0; i < 3; i++) {
        if (point[i] < wall->min_bounds[i] - 1e-6 || point[i] > wall->max_bounds[i] + 1e-6) {
            return 0;
        }
    }
    return 1;
}

// Проекция точки на плоскость стены и проверка границ
int check_wall_collision(const double point[3], const Wall *wall) {
    // Вычисляем расстояние от точки до плоскости
    const double *n = wall->plane.normal;
    double distance = n[0] * point[0] + n[1] * point[1] + n[2] * point[2] - wall->plane.distance;

    // Если точка не на плоскости, пересечения нет
    if (fabs(distance) > 1e-6) return 0;

    // Проверяем, лежит ли проекция точки в пределах стены
    return is_point_in_wall(point, wall);
}

// Основная функция подсчёта коллизий
int collision_counter(const double beg[3], const double end[3], const Box *box, const Room *room) {
    int collisions = 0;
    double dir[3] = {end[0] - beg[0], end[1] - beg[1], end[2] - beg[2]};
    double length = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);

    // Нормализуем направление (для параметризации)
    if (length < 1e-6) return 0; // Отрезок нулевой длины
    dir[0] /= length;
    dir[1] /= length;
    dir[2] /= length;

    // Проверка стен коробки
    for (int i = 0; i < 5; i++) {
        const Wall *wall = &box->walls[i];
        const double *n = wall->plane.normal;
        double d = wall->plane.distance;

        // Параметр t пересечения луча с плоскостью
        double denominator = n[0] * dir[0] + n[1] * dir[1] + n[2] * dir[2];
        if (fabs(denominator) < 1e-6) continue; // Параллельно плоскости

        double t = (d - (n[0] * beg[0] + n[1] * beg[1] + n[2] * beg[2])) / denominator;
        if (t < 0 || t > length) continue; // За пределами отрезка

        // Вычисляем точку пересечения
        double p[3] = {
            beg[0] + dir[0] * t,
            beg[1] + dir[1] * t,
            beg[2] + dir[2] * t
        };

        // Проверяем, лежит ли точка в пределах стены
        if (check_wall_collision(p, wall)) collisions++;
    }

    // Проверка стен комнаты (аналогично)
    for (int i = 0; i < 6; i++) {
        const Wall *wall = &room->walls[i];
        const double *n = wall->plane.normal;
        double d = wall->plane.distance;

        double denominator = n[0] * dir[0] + n[1] * dir[1] + n[2] * dir[2];
        if (fabs(denominator) < 1e-6) continue;

        double t = (d - (n[0] * beg[0] + n[1] * beg[1] + n[2] * beg[2])) / denominator;
        if (t < 0 || t > length) continue;

        double p[3] = {
            beg[0] + dir[0] * t,
            beg[1] + dir[1] * t,
            beg[2] + dir[2] * t
        };

        if (check_wall_collision(p, wall)) collisions++;
    }

    return collisions;
}

//---END_COLLISIONS---///

void copy_matrixes(double v[4][4], const double c[4][4]) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            v[i][j] = c[i][j];
        }
    }
}

void mult_matrixes(double result[4][4], const double a[4][4], const double b[4][4]) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result[i][j] = 0;
            for (int k = 0; k < 4; k++) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void copy_pose(Pose *new_pose, const Pose *pose) {
    new_pose->generalized_coordinates[0] = pose->generalized_coordinates[0];  
    new_pose->generalized_coordinates[1] = pose->generalized_coordinates[1]; 
    new_pose->generalized_coordinates[2] = pose->generalized_coordinates[2]; 
    new_pose->generalized_coordinates[3] = pose->generalized_coordinates[3]; 
    new_pose->generalized_coordinates[4] = pose->generalized_coordinates[4]; 
    new_pose->generalized_coordinates[5] = pose->generalized_coordinates[5]; 
    new_pose->generalized_coordinates[6] = pose->generalized_coordinates[6]; 
}


void initialize_pose(Pose *pose, const Individual *ind) {
    pose->generalized_coordinates[0] = (LINK_0_MIN + LINK_0_MAX) / 2;
    pose->generalized_coordinates[1] = (ind->phenotype[0] + ind->phenotype[1]) / 2;
    pose->generalized_coordinates[2] = (ind->phenotype[2] + ind->phenotype[3]) / 2;
    pose->generalized_coordinates[3] = (ind->phenotype[4] + ind->phenotype[5]) / 2;
    pose->generalized_coordinates[4] = (GRAB_NUT_MAX + GRAB_NUT_MIN) / 2;
    pose->generalized_coordinates[5] = (GRAB_PRC_MAX + GRAB_PRC_MIN) / 2;
    pose->generalized_coordinates[6] = (GRAB_SLF_MAX + GRAB_SLF_MIN) / 2;
}

int dh_parametrs_filling(double **array, const Individual *ind, double *generalized_coordinates) {
	// a_i
	array[0][0] = 0;
	array[1][0] = LINK_1_LEN;
	array[2][0] = ind->phenotype[6];
	array[3][0] = ind->phenotype[7];
	array[4][0] = 0;
	array[5][0] = 0;
	array[6][0] = 0;

	// d_i
	array[0][1] = 0;
	array[1][1] = 0;
	array[2][1] = 0;
	array[3][1] = 0;
	array[4][1] = 0;
	array[5][1] = 0;
	array[6][1] = GRAB_LEN;

	// α_i
	array[0][2] = -M_PI / 2;
	array[1][2] = -M_PI / 2;
	array[2][2] = 0;
	array[3][2] = 0;
	array[4][2] = -M_PI / 2;
	array[5][2] = -M_PI / 2;
	array[6][2] = 0;

	// q_i
	array[0][3] =  0        + generalized_coordinates[0] * M_PI / 180.0;
	array[1][3] = -M_PI / 2 + generalized_coordinates[1] * M_PI / 180.0;
	array[2][3] = -M_PI / 2 + generalized_coordinates[2] * M_PI / 180.0;
	array[3][3] =  0        + generalized_coordinates[3] * M_PI / 180.0;
	array[4][3] =  0        + generalized_coordinates[4] * M_PI / 180.0;
	array[5][3] = -M_PI / 2 + generalized_coordinates[5] * M_PI / 180.0;
	array[6][3] =  0        + generalized_coordinates[6] * M_PI / 180.0;

	return 0;
}

void init_dh_matrix(double A[4][4], double **array, int j) {
   A[0][0] = cos(array[j][3]);
   A[0][1] = -cos(array[j][2]) * sin(array[j][3]);
   A[0][2] = sin(array[j][2]) * sin(array[j][3]);
   A[0][3] = array[j][0] * cos(array[j][3]);

   A[1][0] = sin(array[j][3]);
   A[1][1] = cos(array[j][2]) * cos(array[j][3]);
   A[1][2] = -sin(array[j][2]) * cos(array[j][3]);
   A[1][3] = array[j][0] * sin(array[j][3]);

   A[2][0] = 0;
   A[2][1] = sin(array[j][2]);
   A[2][2] = cos(array[j][2]);
   A[2][3] = array[j][1];

   A[3][0] = 0;
   A[3][1] = 0;
   A[3][2] = 0;
   A[3][3] = 1;
}
