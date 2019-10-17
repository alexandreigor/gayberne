/*
 * SO_specific_func.c
 *
 *  Created on: 09/07/2012
 *      Author: igor
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static inline char* newStringInt(unsigned int size) {
	return (char*) calloc(size, sizeof(char));
}

int copyFile(const char* fileSrc, const char* fileDest) {
	char *comm;
#ifdef OS_WIN
	comm = newStringInt(9 + strlen(fileSrc) + 3 + strlen(fileDest) + 2);
	sprintf(comm, "copy /Y \"%s\" \"%s\"", fileSrc, fileDest);
	//return 1;
#else
	comm = newStringInt(7 + strlen(fileSrc) + 3 + strlen(fileDest) + 2);
	sprintf(comm, "cp -f \'%s\' \'%s\'", fileSrc, fileDest);
#endif
	printf("%s\n", comm);
	int ret = system(comm);
	free(comm);
	return ret;
}

int moveFile(const char* fileSrc, const char* fileDest) {
	char *comm;
#ifdef OS_WIN
	comm = newStringInt(9 + strlen(fileSrc) + 3 + strlen(fileDest) + 2);
	sprintf(comm, "move /Y \'%s\' \'%s\'", fileSrc, fileDest);
#else
	comm = newStringInt(7 + strlen(fileSrc) + 3 + strlen(fileDest) + 2);
	sprintf(comm, "mv -f \'%s\' \'%s\'", fileSrc, fileDest);
#endif
	printf("%s\n", comm);
	int ret = system(comm);
	free(comm);
	return ret;
}
