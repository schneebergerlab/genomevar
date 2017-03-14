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
void parseITX(std::vector<BLOCK> &chromo, char chr[], int num);
void groupITX(std::vector<BLOCK> &chromo, char chr[], int num); //int a, int b);
void writeITX(std::vector<BLOCK> &chromo,int a, int b);

#endif /* GETITX_H_ */
