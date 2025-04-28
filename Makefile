# Компилятор и флаги
CC = gcc
CFLAGS = -Wall -Wextra -pedantic -O2 -g -Iinclude  # -g для отладки, -pedantic для строгости
LDFLAGS = -lm  # Флаги линковщика (математическая библиотека)

# Целевой исполняемый файл
TARGET = genetic_algorithm_a

# Исходные файлы (добавлены simulation.c и config.c)
SRCS = src/EntryPoint.c \
        src/common_functions.c \
        src/logger.c \
        src/simulation.c \
        src/config.c  # Переименованный confg.c → config.c

# Объектные файлы (генерируются в obj/)
OBJS = $(SRCS:src/%.c=obj/%.o)

# Правило по умолчанию
all: $(TARGET)

# Сборка исполняемого файла
$(TARGET): $(OBJS)
	$(CC) -o $@ $(OBJS) $(LDFLAGS)  # LDFLAGS в конце!

# Компиляция .c → .o с созданием obj/ при необходимости
obj/%.o: src/%.c | obj/
	$(CC) $(CFLAGS) -c $< -o $@

# Создание каталога obj/
obj/:
	mkdir -p obj

# Очистка
clean:
	rm -f $(OBJS) $(TARGET)
	rm -rf obj

.PHONY: all clean
