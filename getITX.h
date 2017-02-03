/*
 * getITX.h
 *
 *  Created on: Feb 2, 2017
 *      Author: goel
 */


#ifndef GETITX_H_
#define GETITX_H_
#include "init.h"

//ITX identification
void parseITX(BLOCK *chromo, char chr[], int num);
void groupITX(BLOCK *chromo, char chr[], int num); //int a, int b);
void writeITX(BLOCK *chromo,int a, int b);

#endif /* GETITX_H_ */
