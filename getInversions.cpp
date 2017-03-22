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

	//printf("num:%d\n",num);

	//	printf("INSIDE\n");

	for (int i = 1; i < num-1; i++) {


		int forward[num], backward[num];

		for (int j =1; j<num-1;j++){
			forward[j] =0;
			backward[j] = 0;
		}
		if (strcmp(chromo[i].bchr, chr) == 0 && chromo[i].state != CTX){
			if (chromo[i].dir == -1) {
				/*	 now decide whether the inversion is translocated or not
				 *
				 * Set right B neighbor */

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

						while (chromo[invEndA].dir == 1 || chromo[invEndA].state == CTX ||
								chromo[invEndA].state == SYN || chromo[invEndA].state == ETER ||
								forward[invEndA] == 1 || (chromo[invEndA].bstart > chromo[rightB].bstart && chromo[rightB].state != ETER) ||
								chromo[invEndA].bstart < chromo[leftA].bstart) {
							invEndA--;
						}


						//Run left until next syntenic block
						int leftB = chromo[invEndA].leftBNeighbor;
						if(chromo[leftB].dir == -1) backward[leftB] = 1;

						/** MISSING SOLUTION END OF CHR */
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
						/*for (int j = invBeginA; j <= invEndA; j++) {
							printf("chromo.state: %d\n", chromo[j].state);
						}*/



						if (invBeginA == i) {

							/*	//printf("START\n");


														//printf("Chr: %s \t i: %d\t invBegin: %d\t invEnd: %d\t begin.start: %d \t end.start: %d \t begin.dir: %d \t end.dir: %d\n",
												//				chr,i, invBeginA, invEndA, chromo[invBeginA].bstart, chromo[invEndA].bstart, chromo[invBeginA].dir, chromo[invEndA].dir);

													//	printf("rN: %d\n", chromo[invEndA].rightBNeighbor);




							//printf("Chr: %s \t A.start: %d \n",chr, chromo[i].astart);

							printf("\n%d\t", chromo[i].astart);
							printf("%d\t", chromo[i].aend);
							printf("%d\t", chromo[i].bstart);
							printf("%d\t", chromo[i].bend);
									printf("%d\t", chromo[i].dir);
									printf("%s\t", chromo[i].achr);
									printf("%s\n", chromo[i].bchr);

									printf("invBegin: %d \t invEnd: %d\n",invBeginA, invEndA);


							 */

							fprintf(invOutFile, "#INV ");
							fprintf(invOutFile, "%s %d %d - ", chromo[invBeginA].achr, chromo[invBeginA].astart, chromo[invEndA].aend);
							fprintf(invOutFile, "%s %d %d\n", chromo[invEndA].bchr, chromo[invEndA].bstart, chromo[invBeginA].bend);

							std::vector<int> longestInverted(invEndA - invBeginA +1);
							int invCount = 0 ;


				//			printf("invCount: %d\n",invCount);

				//			printf("invBeginA: %d \t invEndA: %d \n",invBeginA, invEndA);

					//		printf("rightb: %d \t leftA : %d \t invBeginA: %d \t invEndA: %d\n", rightB,leftA, invBeginA, invEndA);




							for (int j = invBeginA; j <= invEndA; j++) {

								longestInverted[j - invBeginA] = 0;


								if (chromo[j].state == SYN) {


									chromo[j].state = SYN_IN_INV;
								//	printf("found SYN \n");
								}else {



									if (chromo[j].state == CTX || chromo[j].state == ITX){
										//printf("found ITX \n");
										continue;
									}
									else {


										if (chromo[j].dir == -1) {


/*

																printf("j.start: %d \t j.end: %d \n"
													"begin.astart : %d begin.aend: %d \t end.astart: %d \t end.aend: %d\n"
													"begin.bstart : %d begin.bend: %d \t end.bstart: %d \t end.bend: %d\n",
													chromo[j].bstart,chromo[j].bend,
													chromo[invBeginA].astart, chromo[invBeginA].aend, chromo[invEndA].astart, chromo[invEndA].aend,
													chromo[invBeginA].bstart, chromo[invBeginA].bend, chromo[invEndA].bstart, chromo[invEndA].bend);*/

									if((chromo[j].bstart > chromo[rightB].bstart || chromo[j].bstart < chromo[leftB].bstart) && chromo[rightB].state != ETER){
											//printf("*****************found ITX_IN *****************************\n");


											//if(chromo[j].bstart > chromo[invBeginA].bstart  ||  chromo[j].bstart < chromo[invEndA].bstart){
									//			printf("*****************found INV_ITX *****************************\n");
												chromo[j].state = INV_ITX;
											}
											else{

									//			printf("found INV \n");
												longestInverted[j-invBeginA] = 1;
												invCount++;
											}
										}
									}
								}

							}



					//		std::cout<<"LOL\n";



						//	printf("invCount: %d\n",invCount);

							//inverted = (BLOCK *)calloc(invCount, sizeof(BLOCK));
							invCount=0;

							for (int j = invBeginA; j <= invEndA; j++){
								if(longestInverted[j-invBeginA] == 1){
									inverted.push_back(chromo[j]);
									invCount++;
								}
							}
						//	std::cout<<"size: "<<inverted.size()<<"   invCount: "<<invCount<<"\n";

					//		std::cout<<chr<<"\t"<<invCount<<"\n";

					//		std::cout<<"inverted.size: "<<inverted.size()<<"\n";

							SYNPATH synPath;

							synPath = invSYN(chr, invCount);


							for(int j =0; j< synPath.maxWeightPathLength; j++){
								inverted[synPath.maxWeightPath[j-1]].state = INV;
								writeBlock(invOutFile, inverted[synPath.maxWeightPath[j-1]]);
							}


							/*for (int j=0; j <inverted.size();++j){
								std::cout<<"state: "<<synPath.invOutFile.state<<"\n";
							}*/


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
							//	free(longestInverted);
							//longestInverted =NULL;
							//	free(inverted);
							//inverted = NULL;


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

	//std::cout<<"length: "<<length<<"\n";

	std::vector<BLOCK> invChrom(length); //= (BLOCK *) calloc(length, sizeof(BLOCK));




	invChrom[0].state = STER;
	invChrom[length-1].state = ETER;



	for (int j =1; j<length-1;j++){
		invChrom[j] = inverted[j-1];
		//		printf("start: %d\t", invChrom[j].bstart);
		//		printf("start: %d\n", invChrom[j].bend);

		invChrom[j].bstart = sum - invChrom[j].bstart;
		invChrom[j].bend = sum - invChrom[j].bend ;
		invChrom[j].dir = (-1) * invChrom[j].dir;

		invChrom[j].bstart = invChrom[j].bstart + invChrom[j].bend;
		invChrom[j].bend = invChrom[j].bstart - invChrom[j].bend;
		invChrom[j].bstart = invChrom[j].bstart - invChrom[j].bend;
	}



//	printf("*****\n");
	/*for (int j=0; j<length;j++){
		printf("%d\t", invChrom[j].astart);
		printf("%d\t", invChrom[j].aend);
		printf("%d\t", invChrom[j].bstart);
		printf("%d\t", invChrom[j].bend);
		printf("%d\t", invChrom[j].dir);
		printf("%s\t", invChrom[j].achr);
		printf("%s\n", invChrom[j].bchr);
	}
*/
	//std::cout<<"invChrom.size: "<<invChrom.size()<<"\n";

	SYNPATH synPath = parseSYN(invChrom,chr,length);

	//free(invChrom);
	//invChrom = NULL;
	//	printf("CHECK_3\n");

//std::cout<<"maxWeightPathLength: "<<synPath.maxWeightPathLength<<"\n";


	return(synPath);


}
