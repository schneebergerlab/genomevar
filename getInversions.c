/*
 * getInversions.c
 *
 *  Created on: Feb 2, 2017
 *      Author: goel
 */

#include "getInversions.h"
#include "getSyn.h"
//#include "SynSearch1.h"



void writeInversions(char chr[], int num, BLOCK *chromo) {

	for (int i = 1; i < num-1; i++) {
		if (strcmp(chromo[i].bchr, chr) == 0 && chromo[i].state != CTX){
			if (chromo[i].dir == -1) {
				/*	 now decide whether the inversion is translocated or not
				 *
				 * Set right B neighbor */

				int rightB = chromo[i].rightBNeighbor;
				while (chromo[rightB].state != SYN && chromo[rightB].state != ETER){
					rightB = chromo[rightB].rightBNeighbor;
				}

				// Set left A neighbor
				int leftA = i - 1;
				while (chromo[leftA].state != STER && chromo[leftA].state != SYN) {
					leftA--;
				}

				if (!(chromo[leftA].state == STER && chromo[rightB].state == ETER)) { //Special case: whole chr inverted, means no inversion as then the alignment dir wrong
					/*if (rightB == -1 || leftA == -1) {
						printf("WARNING: Inv at the end of chr: not solved. Will be reported as ITX.\n");
						if (rightB == -1){
							int leftBlock = i-1;
							int startA = i;
							int startB = blocks[leftBlock].rightBNeighbor;
							int j =0;
							for( j = startB; j > startA; ){

								if(blocks[j].rightBNeighbor == j-1){
									j--;
								}
								else{
									break;
								}
							}

							if (j == i) {
								//printf("INVERSION! invBeginBlock %d invEndBlock %d\n", invBeginA, invEndA);
								fprintf(invOutFile, "#INV ");
								fprintf(invOutFile, "%s %d %d - ", blocks[startA].achr, blocks[startA].astart, blocks[startB].aend);
								fprintf(invOutFile, "%s %d %d\n", blocks[startB].bchr, blocks[startB].bstart, blocks[startA].bend);

								for (int k = startA; k <= startB; k++) {
									if (blocks[k].state == SYN) {
										blocks[k].state = SYN_IN_INV;
									}
									else {
										if (blocks[k].state == CTX) {
											printf("******ANOMALY*********");
											printf("%s %d",chr, blocks[k].astart);}
										else {
											if (blocks[k].dir == -1) {
												blocks[k].state = INV;
												writeBlock(invOutFile, blocks[k]);
											}
										}
									}
								}
								i = startB +1;
							}
						}
						else if(leftA == -1){
							int rightBlock = blocks[i].rightBNeighbor;
							int endA = rightBlock-1;
							int endB =i;

							int startA = i;
							int startB = blocks[rightBlock].rightBNeighbor;
							int j =0;
							for( j = startB; j > startA; ){

								if(blocks[j].rightBNeighbor == j-1){
									j--;
								}
								else{
									break;
								}
							}

							if (j == i) {
								//printf("INVERSION! invBeginBlock %d invEndBlock %d\n", invBeginA, invEndA);
								fprintf(invOutFile, "#INV ");
								fprintf(invOutFile, "%s %d %d - ", blocks[startA].achr, blocks[startA].astart, blocks[startB].aend);
								fprintf(invOutFile, "%s %d %d\n", blocks[startB].bchr, blocks[startB].bstart, blocks[startA].bend);

								for (int k = startA; k <= startB; k++) {
									if (blocks[k].state == SYN) {
										blocks[k].state = SYN_IN_INV;
									}
									else {
										if (blocks[k].state == CTX) {
											printf("******ANOMALY*********");
											printf("%s %d",chr, blocks[k].astart);}
										else {
											if (blocks[k].dir == -1) {
												blocks[k].state = INV;
												writeBlock(invOutFile, blocks[k]);
											}
										}
									}
								}
								i = startB +1;
							}
						}
					}
					else {*/

					if (rightB > i) { // otherwise running back is senseless and no inversion is found

						//Run back on A genome until an inversion is found
						int invEndA = rightB;
						/** if (rightB == -1) { // inv at the end of the chromosome
								invEndA = BLOCK_NUM-1;
							}*/

						while (chromo[invEndA].dir == 1 || chromo[invEndA].state == CTX || chromo[invEndA].state == SYN || chromo[invEndA].state == ETER) {
							invEndA--;
						}

						//Run left until next syntenic block
						int leftB = chromo[invEndA].leftBNeighbor;
						/** MISSING SOLUTION END OF CHR */
						while (chromo[leftB].state != STER && chromo[leftB].state != SYN) {
							leftB = chromo[leftB].leftBNeighbor;
						}

						int invBeginA = leftB;

						while (chromo[invBeginA].dir != -1 || chromo[invBeginA].state == CTX || chromo[invBeginA].state == SYN || chromo[invBeginA].state == STER) {
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
							fprintf(invOutFile, "%s %d %d - ", chromo[invBeginA].achr, chromo[invBeginA].astart, chromo[invEndA].aend);
							fprintf(invOutFile, "%s %d %d\n", chromo[invEndA].bchr, chromo[invEndA].bstart, chromo[invBeginA].bend);

							longestInverted = (int *)calloc((invEndA - invBeginA +1), sizeof(int));
							int invCount=0;

							for (int j = invBeginA; j <= invEndA; j++) {

								longestInverted[j] = 0;
								if (chromo[j].state == SYN) {
									chromo[j].state = SYN_IN_INV;
								}
								else {
									if (chromo[j].state == CTX || chromo[j].state == ITX){
										continue;
									}else {
										if (chromo[j].dir == -1) {
											if(chromo[j].bstart > chromo[invBeginA].bstart  ||  chromo[j].bstart < chromo[invEndA].bstart){
												chromo[j].state = ITX_IN;
											}
											else{
												longestInverted[j] = 1;
												invCount++;
											}
										}
									}
								}
							}
							inverted = (BLOCK *)calloc(invCount, sizeof(BLOCK));
							invCount=0;

							for (int j = invBeginA; j <= invEndA; j++){
								if(longestInverted[j] == 1){
									inverted[invCount] = chromo[j];
									invCount++;
								}
							}

							for(int j =0; j < invCount; j++){
								printf("inv.start: %d\n",inverted[j].astart);
							}


							i = invEndA;
}}}}}}
	//printf("****Finished Writing Inversions****\n");
}


