/*
 * getDuplicates.h
 *
 *  Created on: Mar 16, 2017
 *      Author: goel
 */
#include "init.h"
#ifndef GETDUPLICATES_H_
#define GETDUPLICATES_H_

void parseDUP(std::vector<BLOCK> &blocks, std::vector<BLOCK> &mBlocks);
void filterBlocks(std::vector<BLOCK> &blocks);
void filterBlocks(std::vector<BLOCK> &mBlocks, std::vector<BLOCK> const &blocks);



#endif /* GETDUPLICATES_H_ */
