/*
 * SynSearch.h
 *
 *  Created on: May 1, 2015
 *      Author: schneeberger
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef SYNSEARCH_H_
#define SYNSEARCH_H_

#define MAX_EDGE_NUM 1000
#define MAX_BLOCK_NUM 1000000

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


int BLOCK_NUM = 0;
int CHROMOSOME_NUM = 0;
char CHROMOSOME[4096][4096];

char inputFileName[249];
FILE *inputFile;
FILE *synOutFile;
FILE *ctxOutFile;
FILE *invOutFile;
FILE *itxOutFile;

int maxWeight;
int maxWeightBlock;
int maxWeightPath[MAX_BLOCK_NUM];
int maxWeightPathLength;

//Syntenic Path identification
void setPathWeights(char chr[]);
void setEdges(char chr[]);
int testSynteny(int i, int j);
void backtraceSynPath();
void printSynPath();

//General Helpers
void setEdgesBGenome(char chr[]);
void writeBlock(FILE *file, BLOCK block);
void init(int argc, char *argv[]);
void readInputFile();
void displayHelp(int exitCode);
void displayVersion();

//Inversion identification
void writeInversions(char chr[]);

//CTX identification
void parseCTX();
void writeCTX(int a, int b);

//ITX identification
void parseITX();
void writeITX(int a, int b);

//Debug
void printBlocks();
void printBlock();

const int SYN = 1;
const int SYN_IN_INV = 2;
const int CTX = 3;
const int INV = 4;
const int ITX = 5;



#endif /* SYNSEARCH_H_ */
