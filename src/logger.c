#include "logger.h"

FILE *logFile = NULL;
char logBuffer[4096];


void initLog(const char *filename) {
	logFile = fopen(filename, "w"); // w - clean+write | a - write 
	if (!logFile) {
		perror("Failed to open log file");
		return;
	}
	setvbuf(logFile, logBuffer, _IOFBF, sizeof(logBuffer)); // full bufferisation
	logMessage("TIMESTAMP,EVENT_TYPE,DATA");
}

void closeLog() {
	if (logFile) {
		fclose(logFile);
	}
}


void logMessage(const char *message) {
	if (logFile) {
		fprintf(logFile, "%s\n", message);
	}
}
