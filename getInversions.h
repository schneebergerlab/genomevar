/*
 * getInversions.h
 *
 *  Created on: Feb 2, 2017
 *      Author: goel
 */

#ifndef GETINVERSIONS_H_
#define GETINVERSIONS_H_

#include "init.h"
//Inversion identification

void parseINV(std::vector<BLOCK> &chromo, char chr[], int num);
SYNPATH invSYN(char chr[], int length);

#endif /* GETINVERSIONS_H_ */
