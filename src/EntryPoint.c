
//---HEADERS---//

#include "common_functions.h"
#include "logger.h"
#include "config.h"
#include "structures.h"

//---LIBRARIES---//

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

int main() {
    // LOGS //
    initLog("log.txt");
    logMessage("Program started\n");
    
    // RANDOM REFRESH //
    srand(time(NULL));
    
    // TIME //
    clock_t start, end;
    double total_time;
    start = clock();

    // DEFINITIONS //
    Population pop, new_pop; // Объявление старой и новой популяций
    initialize_population(&pop); // Инициализация текущей популяции  
    // Инициализация окружения 
    Box box;
    initialize_box(&box);
    Room room;
    initialize_room(&room);
    Target target;
    initialize_target(&target);

    // MAIN CYCLE //
    for (int generation = 0; generation < MAX_GENERATIONS; generation++) { // Запускаем поколения
        pop.generation = generation; // Write in generation id
        evaluate_population(&pop, &box, &room, &target, generation); // Оцениваем приспособленность текущей популяции ( каждой особи приписываем показатель приспособленности ) ЗДЕСЬ ПРОИСХОДИТ СИМУЛЯЦИЯ

        for (int i = 0; i < POPULATION_SIZE; i++) { // Заполняем новую популяцию КАЖДЫЙ РАЗ ВЫБИРАЯ НОВЫХ РОДИТЕЛЕЙ И СОЗДАВАЯ НОВОГО НАСЛЕДНИКА
            Individual parent1, parent2;                         // Инициализируем родителей
	    select_parent(&pop, &parent1, &parent2);             // ПровоЖц:дим отбор методом турнирной селекции (заполняем родителей данными)
            Individual child;                                    // Объявляем наследника
            crossover(&parent1, &parent2, &child);               // Скрещиваем (заполняем наследника)
            
	    float k = 0.7f;
	    if (((parent1.fitness + parent2.fitness) / 2) > 200) {
		    k = 0.3f;
            }	    

	    if (i > k * POPULATION_SIZE) {
	        mutate(&child);                                  // Случайно мутируем наследника
	    }

	    
            // mutate(&child); 
            new_pop.individuals[i] = child;                      // Заполняем в i-ю ячейку новой популяции
        }
        replace_population(&pop, &new_pop);
    }
    // TIME STOP //
    end = clock();
    total_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Total time was %.5f seconds or %.5f minutes.\n", total_time, total_time / 60);

    // CLOSING LOGS //
    logMessage("Processing data...");
    closeLog();

    return 0;
}
