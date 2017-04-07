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
#include <string>

#ifndef INIT_H_
#define INIT_H_

#define MAX_EDGE_NUM 10000
#define MAX_BLOCK_NUM 100000

typedef struct block {
	int indexA;
	int indexB;

	int leftBNeighbor; //order on B genome
	int rightBNeighbor;

	std::string achr;
	int astart;
	int aend;
	int alen;

	std::string bchr;
	int bstart;
	int bend;
	int blen;

	int astate;
	int bstate;

//
//	int aDupOf;
//	std::vector<int> aDuplicates;
//
//	int bDupOf;
//	std::vector<int> bDuplicates;

	int dupOf;
	std::vector<int> duplicates;

	float iden;

	int dir;
	int state;

	//A genome (edges describe all potential best syntenic paths)
	int inEdgeNum;
	std::vector<int> inEdge;
	int maxInEdge; // point to the one inEdge that is on the most heavy path


	int outEdgeNum;
	std::vector<int> outEdge;
	//std::map<int,int> outEdge;


	int flagBadOut;
	int weight;
} BLOCK;

//extern BLOCK blocks[MAX_BLOCK_NUM];
extern std::vector<BLOCK> blocks;
extern std::vector<BLOCK> mblocks;

typedef struct synPath{
	int *maxWeightPath;
	int maxWeightPathLength;
}SYNPATH;

typedef struct filteredData{
    std::vector<BLOCK> uniBlocks;
    std::vector<BLOCK> dupBlocks;
} FILTEREDDATA;

extern int BLOCK_NUM;
extern int mBLOCK_NUM;

extern int CHROMOSOME_NUM;
extern int mCHROMOSOME_NUM;

extern char CHROMOSOME[4096][256];
extern char mCHROMOSOME[4096][256];

extern char inputFileName[256];
extern char minputFileName[256];
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
void readInputFile(char *fileName, std::vector<BLOCK> &blocks, int &BLOCK_NUM, int &CHROMOSOME_NUM, char (&CHROMOSOME)[4096][256]);
void displayHelp(int exitCode);
void displayVersion();
void defineUniqueBlocks(std::vector<BLOCK> &blocks, int const &BLOCK_NUM);


//Debug
void printBlocks();
void printBlock();


extern const int RED;				//REDUNDANT
extern const int NA;				//NOT-ASSIGNED
extern const int SYN;				//SUNTENIC
extern const int SYN_IN_INV;		//SYNTENIC INSIDE INVERSION
extern const int CTX;				//CROSS CHROMOSOMAL TRANSLOCATION
extern const int INV;				//INVERSION
extern const int ITX;				//INTER CHROMOSOMAL TRANSLOCATION
extern const int INV_ITX;			//INVERTED INTER CHROMOSOMAL TRANLSOCATION
extern const int ITX_IN_INV;		//INTER CHROMOSOMAL TRANSLOCATION INSIDE INVERSION
extern const int DUP;				//DUPLICATED
extern const int UNI;				//UNIQUE

extern const int TEMP_DUP;			//TEMPORARY DUPLICATE
extern const int STER;				//START TERMINAL
extern const int ETER;				//END TERMINAL


#endif /* INIT_H_ */
