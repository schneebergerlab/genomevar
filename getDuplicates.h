/*
 * getDuplicates.h
 *
 *  Created on: Mar 16, 2017
 *      Author: goel
 */
#include "init.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <fstream>
#ifndef GETDUPLICATES_H_
#define GETDUPLICATES_H_

void parseDUP(std::vector<BLOCK> &blocks, std::vector<BLOCK> &mBlocks, int threshold, FILTEREDDATA &fData);
void filterBlocks(std::vector<BLOCK> &blocks, int threshold);
void filterBlocks(std::vector<BLOCK> &mBlocks, std::vector<BLOCK> const &blocks, int threshold);
void annotateDup(std::vector<BLOCK> &dupBlocks, std::vector<BLOCK> &uniBlocks, int threshold);
void printDup(std::vector<BLOCK> const &dupBlocks, std::vector<BLOCK> const &uniBlocks);
void filtermUni(std::vector<BLOCK> &blocks, int threshold);
std::vector<int> getStateCount(std::vector<BLOCK> const &blocks);
void sortBlocks(std::vector<BLOCK>& blocks);

#endif /* GETDUPLICATES_H_ */
