/*
 * getInversions.c
 *
 *  Created on: Feb 2, 2017
 *      Author: goel
 */

#include "getInversions.h"
#include "getSyn.h"
#include <vector>
//#include "SynSearch1.h"

int *longestInverted;
std::vector<BLOCK> inverted;

void parseINV(std::vector<BLOCK> &chromo, char chr[], int num) {

  //  std::cout<<"START\n";
	for (int i = 1; i < num-1; i++) {

		int forward[num], backward[num];
		for (int j =1; j<num-1;j++){
			forward[j] =0;
			backward[j] = 0;
		}

		if (chromo[i].bchr.compare(chr) == 0 && chromo[i].state != CTX){
			if (chromo[i].dir == -1) {
				int rightB = chromo[i].rightBNeighbor;
				if(chromo[rightB].dir == -1) forward[rightB] = 1;
				while (chromo[rightB].state != SYN && chromo[rightB].state != ETER){
					rightB = chromo[rightB].rightBNeighbor;
					if(chromo[rightB].dir == -1) forward[rightB] = 1;
				}
				int leftA = i - 1;
				while (chromo[leftA].state != STER && chromo[leftA].state != SYN) {
					leftA--;
				}
				if (!(chromo[leftA].state == STER && chromo[rightB].state == ETER)) {
					//Special case: whole chr inverted, means no inversion as then the alignment dir wrong
					if (rightB > i) { // otherwise running back is senseless and no inversion is found
						//Run back on A genome until an inversion is found
						int invEndA = rightB;
						while (chromo[invEndA].dir == 1 || chromo[invEndA].state == CTX ||
								chromo[invEndA].state == SYN || chromo[invEndA].state == ETER ||
								forward[invEndA] == 1 || (chromo[invEndA].bstart > chromo[rightB].bstart && chromo[rightB].state != ETER) ||
								chromo[invEndA].bstart < chromo[leftA].bstart) {
							invEndA--;
						}

						//Run left until next syntenic block
						int leftB = chromo[invEndA].leftBNeighbor;
						if(chromo[leftB].dir == -1) backward[leftB] = 1;

						/** MISSING SOLUTION END OF CHR **/
						while (chromo[leftB].state != STER && chromo[leftB].state != SYN) {
							leftB = chromo[leftB].leftBNeighbor;
							if(chromo[leftB].dir == -1) backward[leftB] = 1;

						}

						int invBeginA = leftB;
						while (chromo[invBeginA].dir != -1 || chromo[invBeginA].state == CTX ||
								chromo[invBeginA].state == SYN || chromo[invBeginA].state == STER ||
								backward[invBeginA] == 1 ||  (chromo[invBeginA].bstart > chromo[rightB].bstart && chromo[rightB].state != ETER) ||
								chromo[invBeginA].bstart < chromo[leftA].bstart) {
							invBeginA++;
							if(invBeginA > num-1){
								break;
							}
						}
						// If the end of the circuit equals the beginning: Inversion found
						// This can include rearranged inversions in the inversion.
						// Actually the best would be to take inversions and run SnyPathFinder within them (respecting the inversion)
						// This would give the optimal inversion and report the ITX in the inversion

						if (invBeginA == i) {
							fprintf(invOutFile, "#INV ");
							fprintf(invOutFile, "%s %d %d - ", chromo[invBeginA].achr.c_str(), chromo[invBeginA].astart, chromo[invEndA].aend);
							fprintf(invOutFile, "%s %d %d\n", chromo[invEndA].bchr.c_str(), chromo[invEndA].bstart, chromo[invBeginA].bend);
							std::vector<int> longestInverted(invEndA - invBeginA +1);
							int invCount = 0 ;
							for (int j = invBeginA; j <= invEndA; j++) {
								longestInverted[j - invBeginA] = 0;
								if (chromo[j].state == SYN) {
									chromo[j].state = SYN_IN_INV;
								}else {
									if (chromo[j].state == CTX || chromo[j].state == ITX){
										continue;
									}
									else {
										if (chromo[j].dir == -1) {
                                            if((chromo[j].bstart > chromo[rightB].bstart || chromo[j].bstart < chromo[leftB].bstart) && chromo[rightB].state != ETER){
                                            //if(chromo[j].bstart > chromo[invBeginA].bstart  ||  chromo[j].bstart < chromo[invEndA].bstart){
                                                chromo[j].state = INV_ITX;
                                            }
											else{
												longestInverted[j-invBeginA] = 1;
												invCount++;
											}
										}
									}
								}
							}
							invCount=0;
							for (int j = invBeginA; j <= invEndA; j++){
								if(longestInverted[j-invBeginA] == 1){
									inverted.push_back(chromo[j]);
									invCount++;
								}
							}
							SYNPATH synPath;
							synPath = invSYN(chr, invCount);
							for(int j =0; j< synPath.maxWeightPathLength; j++){
								inverted[synPath.maxWeightPath[j-1]].state = INV;
								writeBlock(invOutFile, inverted[synPath.maxWeightPath[j-1]]);
							}
							invCount =0;
							for (int j = invBeginA; j <= invEndA; j++) {
								if(longestInverted[j-invBeginA] == 1){
									if(inverted[invCount].state == INV){
										chromo[j].state = INV;
									}
									else{
										chromo[j].state = ITX_IN_INV;
									}
									invCount++;
								}
							}
							i = invEndA;
							inverted.clear();
							free(synPath.maxWeightPath);
                        }
                    }
                }
            }
        }
    }
}

SYNPATH invSYN(char chr[], int length){
	int sum = inverted[0].bend + inverted[length-1].bstart;
	length = length+2;
	std::vector<BLOCK> invChrom(length);
	invChrom[0].state = STER;
	invChrom[length-1].state = ETER;
	for (int j =1; j<length-1;j++){
		invChrom[j] = inverted[j-1];
		invChrom[j].bstart = sum - invChrom[j].bstart;
		invChrom[j].bend = sum - invChrom[j].bend ;
		invChrom[j].dir = (-1) * invChrom[j].dir;

		invChrom[j].bstart = invChrom[j].bstart + invChrom[j].bend;
		invChrom[j].bend = invChrom[j].bstart - invChrom[j].bend;
		invChrom[j].bstart = invChrom[j].bstart - invChrom[j].bend;
	}
	/*for (int j=0; j<length;j++){
		printf("%d\t", invChrom[j].astart);
		printf("%d\t", invChrom[j].aend);
		printf("%d\t", invChrom[j].bstart);
		printf("%d\t", invChrom[j].bend);
		printf("%d\t", invChrom[j].dir);
		printf("%s\t", invChrom[j].achr);
		printf("%s\n", invChrom[j].bchr);
	}*/
	SYNPATH synPath = parseSYN(invChrom,chr,length);
	return(synPath);
}
