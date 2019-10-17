/*
 * SO_specific_func.h
 *
 *  Created on: 09/07/2012
 *      Author: igor
 */

#ifndef SO_SPECIFIC_FUNC_H_
#define SO_SPECIFIC_FUNC_H_

#if defined(_WIN32) || defined(WIN32) || defined(__CYGWIN__) || defined(__MINGW32__) || defined(__BORLANDC__)
#define OS_WIN
#endif

int copyFile(const char* fileSrc, const char* fileDest);
int moveFile(const char* fileSrc, const char* fileDest);

#endif /* SO_SPECIFIC_FUNC_H_ */
