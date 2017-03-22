/*
 * SynSearch1.c
 *
 *  Created on: Feb 2, 2017
 *      Author: goel
 */

/*
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
*/
#include "init.h"
#include<iostream>

//using namespace std;

int BLOCK_NUM=0;
int mBLOCK_NUM=0;

int CHROMOSOME_NUM=0;
int mCHROMOSOME_NUM=0;

std::vector<BLOCK> blocks;
std::vector<BLOCK> mblocks;

char CHROMOSOME[4096][256];
char mCHROMOSOME[4096][256];

char inputFileName[256];
char minputFileName[256];
FILE *inputFile;
FILE *synOutFile;
FILE *ctxOutFile;
FILE *invOutFile;
FILE *itxOutFile;
FILE *dupOutFile;

const int RED = -1;
const int NA = 0;
const int SYN = 1;
const int SYN_IN_INV = 2;
const int CTX = 3;
const int INV = 4;
const int ITX = 5;
const int INV_ITX = 6;
const int ITX_IN_INV =7;
const int DUP = 8;
const int UNI = 9;

const int TEMP_DUP = 100;
const int STER = 101;
const int ETER = 102;

//BLOCK *mBlocks = (BLOCK *) calloc(100, sizeof(BLOCK));

void init(int argc, char *argv[]) {
	char vers[] = "--version";

	if (argc == 1) {
		displayHelp(EXIT_SUCCESS);
	}
	if ((argc == 2) && (argv[1][0]=='-')) {
		if (strcmp(argv[1],&vers[0])==0) {
			displayVersion();
			exit(EXIT_SUCCESS);
		}
		displayHelp(1);
	}
	if(argc == 3 and strcmp(argv[1], "-1")==0){
		//call_1(int argc, char *argv[]);
		strcpy(inputFileName, argv[2]);
		readInputFile(inputFileName, blocks, BLOCK_NUM, CHROMOSOME_NUM, CHROMOSOME);
		//filterBlocks(blocks, BLOCK_NUM);

	}
	if(argc == 5 and strcmp(argv[1], "-1")==0 and strcmp(argv[3],"-m")==0){
		//call_m(int argc, char *argv[]);
		strcpy(inputFileName, argv[2]);
		strcpy(minputFileName, argv[4]);
		readInputFile(inputFileName, blocks, BLOCK_NUM, CHROMOSOME_NUM, CHROMOSOME);
		readInputFile(minputFileName, mblocks, mBLOCK_NUM, mCHROMOSOME_NUM, mCHROMOSOME);
		std::cout<<"size: "<<BLOCK_NUM<<" chromo_num: "<<CHROMOSOME_NUM<<"\n";
		std::cout<<"size: "<<mBLOCK_NUM<<" chromo_num: "<<CHROMOSOME_NUM<<"\n";
	//	filterBlocks(blocks, BLOCK_NUM);
	//	filterBlocks(mblocks, mBLOCK_NUM);
	}


	// open output files:
	char outfilename[] = "SynSearch.syn.txt";
	if ((synOutFile = fopen(outfilename, "w")) == NULL) {
		printf("Cannot open output file\n");
		exit(1);
	}

	char outfilename2[] = "SynSearch.ctx.txt";
	if ((ctxOutFile = fopen(outfilename2, "w")) == NULL) {
		printf("Cannot open output file\n");
		exit(1);
	}

	char outfilename3[] = "SynSearch.inv.txt";
	if ((invOutFile = fopen(outfilename3, "w")) == NULL) {
		printf("Cannot open output file\n");
		exit(1);
	}

	char outfilename4[] = "SynSearch.itx.txt";
	if ((itxOutFile = fopen(outfilename4, "w")) == NULL) {
		printf("Cannot open output file\n");
		exit(1);
	}

	char outfilename5[] = "SynSearch.dup.txt";
	if ((dupOutFile = fopen(outfilename5, "w")) == NULL) {
		printf("Cannot open output file\n");
		exit(1);
	}
std::cout<<"HERE\n";
freopen("log.txt", "w", stderr);
//	if ((log = std::freopen("log.txt", "w", stderr)) == NULL) {
//		printf("Cannot open output file\n");
//		exit(1);
//	}
	std::cout<<"EXITING INIT\n";


}

void displayHelp(int exitCode) {

	printf("USAGE: SynSearch [-options] <coords-file>\n");
	printf("\n");
	printf("options:\n");
	printf("    --version       displays the version number and exits\n");
	printf("\n");
	exit(exitCode);

}

void displayVersion() {
	//printf("SynSearch, version 0.01, 01.05.2015\n");
	//printf("SynSearch, version 0.02, 09.02.2017\n");
	printf("SynSearch, version 0.03, 31.01,2017\n");
}

void readInputFile(char *fileName, std::vector<BLOCK> &blocks, int &BLOCK_NUM, int &CHROMOSOME_NUM, char (&CHROMOSOME)[4096][256]) {

	if ((inputFile = fopen(fileName, "r")) == NULL) {
		printf("Cannot open input file\n");
		exit(1);
	}

	std::cout<<"READING\n";

	int i1, i2, i3, i4, i5, i6, i8, i9;
	float f7;
	char c10[256], c11[256];
	int ret=0, adir, bdir;
	int c = 0;
	int line = 0;
	char currchr[4096];
	int currpos = 0;

	//read in file


	while((ret = fscanf(inputFile, "%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%s\t%s\n",
			&i1, &i2, &i3, &i4, &i5, &i6, &f7, &i8, &i9, &c10[0], &c11[0])) == 11) {

		line++;

		BLOCK new_block;

		if(i1<i2){
			new_block.astart = i1;
			new_block.aend = i2;
			adir = 1;
		}
		else{
			printf("Error in input file. Blocks in the first genome need to be in fwd direction (line %d)\n", line);
			exit(1);
		}

		if (i3 < i4) {
			new_block.bstart = i3;
			new_block.bend = i4;
			bdir = 1;
		} else {
			new_block.bstart = i4;
			new_block.bend = i3;
			bdir = -1;
		}


		new_block.alen = new_block.aend - new_block.astart + 1;
		new_block.blen = new_block.bend - new_block.bstart + 1;

		new_block.iden = f7;

		if (adir == bdir) {
			new_block.dir = 1;
		} else {
			new_block.dir = -1;
		}
		strcpy(new_block.achr, c10);
		strcpy(new_block.bchr, c11);

		new_block.astate = NA;
		new_block.bstate = NA;
		new_block.state = NA;

		c++;

		if (strcmp(c10, currchr) == 0) { // same chr
			if (i1 < currpos) {
				printf("Error in input file. Positions in the first genome need to be sorted and larger 0 (line %d)\n", line);
				exit(1);
			}
			currpos = i1;
		}
		else { // different chromosome
			for (int i = 0; i < CHROMOSOME_NUM; i++) { // there already?
				if (strcmp(CHROMOSOME[i], c10) == 0) {
					printf("Error in input file. Chromosomes are not sorted (line %d).\n", line);
					exit(1);
				}
			}
			strcpy(CHROMOSOME[CHROMOSOME_NUM], c10);
			CHROMOSOME_NUM++;
			strcpy(currchr, c10);
			currpos = i1;
		}
		blocks.push_back(new_block);
	}
	BLOCK_NUM = c;

	if (ret != EOF) {
		printf("Input file error\n");
		exit(1);
	}
	printf("FINSHED READING\n");
	fclose(inputFile);
}

void writeBlock(FILE *file, BLOCK block) {
	fprintf(file, "%d\t", block.astart);
	fprintf(file, "%d\t", block.aend);
	fprintf(file, "%d\t", block.bstart);
	fprintf(file, "%d\t", block.bend);
	fprintf(file, "%d\t", block.dir);
	fprintf(file, "%s\t", block.achr);
	fprintf(file, "%s\n", block.bchr);
}

void setEdgesBGenome(std::vector<BLOCK> &chromo, char chr[], int num) {

	int first = -1, last = -1;

	for (int i = 1; i < num; i++) {
		chromo[i].leftBNeighbor = -1;
		chromo[i].rightBNeighbor = -1;
		if (strcmp(chromo[i].bchr, chr) == 0 && chromo[i].state != CTX){ // only one chr at a time and only on one chr

			// This is not a good sorting algorithm, but as the list is almost sorted this should work.

			//initiate list:
			if (last == -1) {
				first = i;
				last = i;
			}
			// if there is at least one element:
			else {
				int current = last;
				while (current != -1 && chromo[current].bstart > chromo[i].bstart) {
					current = chromo[current].leftBNeighbor;
				}
				if (current == -1) { // i sits at the beginning
					chromo[i].rightBNeighbor = first;
					chromo[first].leftBNeighbor = i;
					first = i;
				}
				else {
					chromo[i].leftBNeighbor = current;
					chromo[i].rightBNeighbor = chromo[current].rightBNeighbor;

					chromo[chromo[i].leftBNeighbor].rightBNeighbor = i;
					if (chromo[i].rightBNeighbor != -1) {
						chromo[chromo[i].rightBNeighbor].leftBNeighbor = i;
					}
					else {
						last = i;
					}
				}
			}
		}
	}
	if(last != -1){
		chromo[0].rightBNeighbor = first;
		chromo[0].leftBNeighbor = -1;
		chromo[first].leftBNeighbor = 0;
		chromo[num-1].leftBNeighbor = last;
		chromo[num-1].rightBNeighbor = -1;
		chromo[last].rightBNeighbor = num-1;
	}
}
