/*
 * getDuplicates.h
 *
 *  Created on: Mar 16, 2017
 *      Author: goel
 */
#include "init.h"
#include <iostream>
#include <algorithm>
#ifndef GETDUPLICATES_H_
#define GETDUPLICATES_H_

FILTEREDDATA parseDUP(std::vector<BLOCK> &blocks, std::vector<BLOCK> &mBlocks, int threshold);
void filterBlocks(std::vector<BLOCK> &blocks, int threshold);
void filterBlocks(std::vector<BLOCK> &mBlocks, std::vector<BLOCK> const &blocks, int threshold);



#endif /* GETDUPLICATES_H_ */
