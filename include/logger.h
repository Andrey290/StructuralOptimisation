#ifndef LOGGER_H
#define LOGGER_H

#include <stdio.h>
#include <time.h>

extern FILE *logFile; // Global var (should be)

void initLog(const char *filename);
void closeLog();
void logMessage(const char *message);

#define LOG_EVENT(event_type, format, ...) { \
    char buffer[512]; \
    time_t now = time(NULL); \
    struct tm *tm_info = localtime(&now); \
    strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", tm_info); \
    char log_line[512]; \
    snprintf(log_line, sizeof(log_line), "%s,%s," format, buffer, event_type, ##__VA_ARGS__); \
    logMessage(log_line); \
}

#endif
