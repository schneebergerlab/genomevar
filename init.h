/*
 * SynSearch.h
 *
 *  Created on: May 1, 2015
 *      Author: schneeberger
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef INIT_H_
#define INIT_H_

#define MAX_EDGE_NUM 1000
#define MAX_BLOCK_NUM 100000

typedef struct block {
	int indexA;
	int indexB;

	int leftBNeighbor; //order on B genome
	int rightBNeighbor;

	char achr[4096];
	int astart;
	int aend;
	int alen;

	char bchr[4096];
	int bstart;
	int bend;
	int blen;

	int dir;
	int state;

	//A genome (edges describe all potential best syntenic paths)
	int inEdgeNum;
	int inEdge[MAX_EDGE_NUM];
	int maxInEdge; // point to the one inEdge that is on the most heavy path
	

	int outEdgeNum;
	int outEdge[MAX_EDGE_NUM];

	char flagBadOut;
	int weight;
} BLOCK;

BLOCK blocks[MAX_BLOCK_NUM];

typedef struct synPath{
	int *maxWeightPath;
	int maxWeightPathLength;
}SYNPATH;

extern int BLOCK_NUM;
extern int CHROMOSOME_NUM;
char CHROMOSOME[4096][4096];

char inputFileName[249];
FILE *inputFile;
FILE *synOutFile;
FILE *ctxOutFile;
FILE *invOutFile;
FILE *itxOutFile;


//General Helpers
//void setEdgesBGenome(char chr[]);
void setEdgesBGenome(BLOCK *chromo, char chr[], int num);
void writeBlock(FILE *file, BLOCK block);
void init(int argc, char *argv[]);
void readInputFile();
void displayHelp(int exitCode);
void displayVersion();


//Debug
void printBlocks();
void printBlock();


extern const int SYN;
extern const int SYN_IN_INV;
extern const int CTX;
extern const int INV;
extern const int ITX;
extern const int INV_ITX;
extern const int ITX_IN_INV;
extern const int STER;
extern const int ETER;


#endif /* INIT_H_ */
