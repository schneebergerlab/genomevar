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
SYNPATH parseSYN(BLOCK *chromo,char chr[], int num);
void setPathWeights(BLOCK *chromo, char chr[], int num);
void setEdges(BLOCK *chromo, char chr[], int num);
void backtraceSynPath(BLOCK *chromo, char chr[], int num);
int testSynteny(int i, int j, BLOCK *chromo);

void printSynPath(BLOCK *chromo, SYNPATH synPath);

#endif /* GETSYN_H_ */
