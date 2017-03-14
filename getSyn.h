/*
 * getSyn.h
 *
 *  Created on: Feb 2, 2017
 *      Author: goel
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<iostream>
#include<ctime>

#include "init.h"
#ifndef GETSYN_H_
#define GETSYN_H_

//Syntenic Path identification
SYNPATH parseSYN(std::vector<BLOCK> &chromo,char chr[], int num);
void setPathWeights(std::vector<BLOCK> &chromo, char chr[], int num);
void setEdges(std::vector<BLOCK> &chromo, char chr[], int num);
void backtraceSynPath(std::vector<BLOCK> &chromo, char chr[], int num);
int testSynteny(int i, int j, std::vector<BLOCK> &chromo);

void printSynPath(std::vector<BLOCK> &chromo, SYNPATH synPath);

#endif /* GETSYN_H_ */
