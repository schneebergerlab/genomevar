/*
 * getSyn.c
 *
 *  Created on: Feb 2, 2017
 *      Author: goel
 */

#include "getSyn.h"

//BLOCK *chromo;
int maxWeight =0;
int maxWeightBlock=0;
int *maxWeightPath;
int maxWeightPathLength=0;


SYNPATH parseSYN(std::vector<BLOCK> &chromo, char chr[], int num){


	printf("********NUM: %d\n", num);
	maxWeightPath = (int *) calloc(num, sizeof(int));
	//test for commit

	setEdgesBGenome(chromo, chr, num);

	// build all possible paths
	setEdges(chromo, chr, num);


	// calculate cummulative weigths for each path
	// while only the heaviest path is followed



	setPathWeights(chromo, chr, num);

	backtraceSynPath(chromo, chr, num);

	//printSynPath(chromo, chr, num);

	SYNPATH synPath;
	synPath.maxWeightPath = maxWeightPath;
	synPath.maxWeightPathLength = maxWeightPathLength;

	//	printf("FINISHED SYN\n");
	maxWeight =0;
	maxWeightBlock=0;
	maxWeightPath = NULL;
	maxWeightPathLength=0;
	return(synPath);
}

void setEdges(std::vector<BLOCK> &chromo, char chr[], int num) {

	for (int i = num-2 ; i > 0; i--) {
		printf("i: %d\n",i);
		if (strcmp(chromo[i].achr, chr) == 0 && chromo[i].dir == 1 && chromo[i].state != CTX){ // only one chromosome at a time

			//printf("i: %d \t Weight: %d \n",i, chromo[i].weight);

			maxWeightBlock = i;
			maxWeight = chromo[i].alen;
			int setEdge = 0;
			int j = i+1;
			chromo[i].flagBadOut = 1;
			while (setEdge == 0 && j < num-1) { // block j needs to existpr
				if (strcmp(chromo[j].achr, chr) == 0 && chromo[j].dir == 1) {

					if (testSynteny(i, j, chromo)) {
						//printf("i: %d \t j: %d\n", i,j);

						chromo[i].outEdge[chromo[i].outEdgeNum] = j;
						chromo[i].outEdgeNum++;
						chromo[j].inEdge[chromo[j].inEdgeNum] = i;
						chromo[j].inEdgeNum++;
						// set badOut to 0 if it is the next block that is linked
						if (i+1 == j && chromo[i].rightBNeighbor == j) {
							chromo[i].flagBadOut = 0;
							setEdge = 1;
						}
						// Only if the linked block is not badOut stop
						// This is relevant if multiple translocations follow each other: multiple outEdges will be set for the block
						/*	if (chromo[i].flagBadOut == 0) {				//Confirm With Korbinian
								printf(" **** i : %d\t j: %d\n",i,j);
								setEdge = 1;
							}
						}}*/
					}
				}
				j++;
			}
		}

	}


}

int testSynteny(int i, int j, std::vector<BLOCK> &chromo) {
	if (chromo[i].astart < chromo[j].astart && chromo[i].bstart < chromo[j].bstart && strcmp(chromo[i].achr, chromo[j].bchr)==0) {
		return 1;
	} else {
		return 0;
	}}

void setPathWeights(std::vector<BLOCK> &chromo, char chr[], int num) {

	for (int i = 1; i < num-1; i++) {
		if (strcmp(chromo[i].achr, chr) == 0 && chromo[i].dir == 1 && chromo[i].state != CTX){ // only one chromosome at a time and only fwd chromo
			chromo[i].weight = chromo[i].alen;
			for (int j = 0; j < chromo[i].inEdgeNum; j++) {
				if (chromo[i].weight < chromo[chromo[i].inEdge[j]].weight + chromo[i].alen) {
					chromo[i].weight = chromo[chromo[i].inEdge[j]].weight + chromo[i].alen;
					chromo[i].maxInEdge = chromo[i].inEdge[j];
					if (chromo[i].weight > maxWeight) {
						printf("i: %d \t j: %d\n", i,j);

						maxWeightBlock = i;
						maxWeight = chromo[i].weight;
					}
				}
			}
		}
	}
}

void backtraceSynPath( std::vector<BLOCK> &chromo, char chr[], int num){
	int maxWeightPathReverse[num-2];
	int length = 0;
	int cblock = maxWeightBlock;

	while (chromo[cblock].inEdgeNum > 0 && strcmp(chromo[cblock].achr, chr)==0) {
		maxWeightPathReverse[length] = cblock;
		cblock = chromo[cblock].maxInEdge;
		length++;
	}

	if(strcmp(chromo[cblock].achr,chr)==0){
		maxWeightPathReverse[length] = cblock;
		length++;
	}

	for (int i = length-1; i >= 0; i--) {
		maxWeightPath[length-i-1] = maxWeightPathReverse[i];
		chromo[maxWeightPath[length-i-1]].state = SYN;
	}

	maxWeightPathLength = length;
}

void printSynPath(std::vector<BLOCK> &chromo, SYNPATH synPath) {

	maxWeightPathLength = synPath.maxWeightPathLength;
	maxWeightPath = synPath.maxWeightPath;

	std::cout<<"maxWeightPathLength: "<<maxWeightPathLength<<"\n";

	int in = 0, s = -1;
	for (int p = 0; p < maxWeightPathLength; p++) {
		int i = maxWeightPath[p];
		if (in == 0) {

			in = 1;
			s = p;
		}
		else {
			// Check on A genome and on B genome if the block before is
			// the neighbor block, then they are in one synentic run.
			if (i-1 != maxWeightPath[p-1] || chromo[i].leftBNeighbor != maxWeightPath[p-1]) { // new region
				fprintf(synOutFile, "#SYN ");
				fprintf(synOutFile, "%s %d %d - ", chromo[maxWeightPath[s]].achr, chromo[maxWeightPath[s]].astart, chromo[maxWeightPath[p-1]].aend);
				fprintf(synOutFile, "%s %d %d\n", chromo[maxWeightPath[s]].bchr, chromo[maxWeightPath[s]].bstart, chromo[maxWeightPath[p-1]].bend);
				for (int j = s; j < p; j++) {
					writeBlock(synOutFile, chromo[maxWeightPath[j]]);
				}
				s = p;
				in = 1;
			}
		}
	}

	if (in == 1) {
		int e = maxWeightPathLength-1;
		fprintf(synOutFile, "#SYN ");
		fprintf(synOutFile, "%s %d %d - ", chromo[maxWeightPath[s]].achr, chromo[maxWeightPath[s]].astart, chromo[maxWeightPath[e]].aend);
		fprintf(synOutFile, "%s %d %d\n", chromo[maxWeightPath[s]].bchr, chromo[maxWeightPath[s]].bstart, chromo[maxWeightPath[e]].bend);
		for (int j = s; j < maxWeightPathLength; j++) {
			writeBlock(synOutFile, chromo[maxWeightPath[j]]);
		}
	}

	maxWeightPath = NULL;
	maxWeightPathLength = 0;
}