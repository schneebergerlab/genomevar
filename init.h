/*
 * SynSearch.h
 *
 *  Created on: May 1, 2015
 *      Author: schneeberger
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#ifndef INIT_H_
#define INIT_H_

#define MAX_EDGE_NUM 10000
#define MAX_BLOCK_NUM 100000

typedef struct block {
	int indexA;
	int indexB;

	int leftBNeighbor; //order on B genome
	int rightBNeighbor;

	char achr[256];
	int astart;
	int aend;
	int alen;

	char bchr[256];
	int bstart;
	int bend;
	int blen;

	float iden;

	int dir;
	int state;

	//A genome (edges describe all potential best syntenic paths)
	int inEdgeNum;
	int inEdge[MAX_EDGE_NUM];
	int maxInEdge; // point to the one inEdge that is on the most heavy path
	

	int outEdgeNum;
	int outEdge[MAX_EDGE_NUM];
	//std::map<int,int> outEdge;


	int flagBadOut;
	int weight;
} BLOCK;

//extern BLOCK blocks[MAX_BLOCK_NUM];
extern std::vector<BLOCK> blocks;
extern std::vector<BLOCK> mBlocks;

typedef struct synPath{
	int *maxWeightPath;
	int maxWeightPathLength;
}SYNPATH;

extern int BLOCK_NUM;
extern int mBLOCK_NUM;

extern int CHROMOSOME_NUM;
extern int mCHROMOSOME_NUM;

extern char CHROMOSOME[4096][4096];
extern char mCHROMOSOME[4096][4096];

extern char inputFileName[249];
extern char minputFileName[249];
extern FILE *inputFile;
extern FILE *synOutFile;
extern FILE *ctxOutFile;
extern FILE *invOutFile;
extern FILE *itxOutFile;
extern FILE *dupOutFile;


//General Helpers
//void setEdgesBGenome(char chr[]);
void setEdgesBGenome(std::vector<BLOCK> &chromo, char chr[], int num);
void writeBlock(FILE *file, BLOCK block);
void init(int argc, char *argv[]);
void readInputFile(char *fileName, std::vector<BLOCK> &blocks);
void displayHelp(int exitCode);
void displayVersion();
void filterBlocks();


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
extern const int DUP;
extern const int TEMP_DUP;
extern const int STER;
extern const int ETER;


#endif /* INIT_H_ */
