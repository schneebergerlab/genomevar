/*
 * getSyn.h
 *
 *  Created on: Feb 2, 2017
 *      Author: goel
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "init.h"
#ifndef GETSYN_H_
#define GETSYN_H_

//Syntenic Path identification
void parseSYN(BLOCK *chromo,char chr[], int num);
void setPathWeights(char chr[], int num);
void setEdges(char chr[], int num);
void backtraceSynPath(char chr[], int num);
void printSynPath();
int testSynteny(int i, int j, BLOCK *chromo);

#endif /* GETSYN_H_ */
