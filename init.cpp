/*
 * SynSearch1.c
 *
 *  Created on: Feb 2, 2017
 *      Author: goel
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "init.h"
#include<iostream>

//using namespace std;

int BLOCK_NUM=0;
int CHROMOSOME_NUM=0;

std::vector<BLOCK> blocks;

char CHROMOSOME[4096][4096];
char mCHROMOSOME[4096][4096];

 char inputFileName[249];
 char minputFileName[249];
 FILE *inputFile;
 FILE *synOutFile;
 FILE *ctxOutFile;
 FILE *invOutFile;
 FILE *itxOutFile;
 FILE *dupOutFile;

const int SYN = 1;
const int SYN_IN_INV = 2;
const int CTX = 3;
const int INV = 4;
const int ITX = 5;
const int INV_ITX = 6;
const int ITX_IN_INV =7;
const int DUP = 8;

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

  strcpy(inputFileName, argv[argc-1]);
 // strcpy(minputFileName, argv[argc-1]);

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


	readInputFile(inputFileName, blocks);

	std::cout<<"size: "<<blocks.size()<<"\n";

	filterBlocks();

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


//std::vector<BLOCK> readInputFile(char *fileName, std::vector<BLOCK> blocks) {
void readInputFile(char *fileName, std::vector<BLOCK> &blocks) {
  if ((inputFile = fopen(fileName, "r")) == NULL) {
	  printf("Cannot open input file\n");
    exit(1);
  }

  int i1, i2, i3, i4, i5, i6, i8, i9;
  float f7;
  char c10[4096], c11[4096];
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

	 	  new_block.iden = f7;

	 	  if (adir == bdir) {
	 		  new_block.dir = 1;
	 	  } else {
	 		  new_block.dir = -1;
	 	  }

	 	  strcpy(new_block.achr, c10);
	 	 strcpy(new_block.bchr, c11);

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
/*
	  //read in block from A genome
	  if (i1 < i2) {
		  blocks[c].astart = i1;
		  blocks[c].aend = i2;
		  adir = 1;
      	  } else {
		  printf("Error in input file. Blocks in the first genome need to be in fwd direction (line %d)\n", line);
		  exit(1);
      	  }


	  // read in block from B genome: consider direction
	  if (i3 < i4) {
          blocks[c].bstart = i3;
          blocks[c].bend = i4;
          bdir = 1;
      } else {
          blocks[c].bstart = i4;
          blocks[c].bend = i3;
          bdir = -1;
      }

	  blocks[c].alen = blocks[c].aend - blocks[c].astart + 1;

	  blocks[c].iden = f7;

	  if (adir == bdir) {
		  blocks[c].dir = 1;
	  } else {
		  blocks[c].dir = -1;
	  }

	  //record chromosome/scaffold/contig of that block in both genomes
	  strcpy(blocks[c].achr, c10);
	  strcpy(blocks[c].bchr, c11);

	  // increment block counter
	  c++;

	  // collect chromosome identifier:
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

	  */
  }
  BLOCK_NUM = c;

  if (ret != EOF) {
      printf("Input file error\n");
      exit(1);
  }
printf("FINSHED READING\n");
  fclose(inputFile);

 // return(blocks);

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

	//printf("CHROMOSOME: %s\n",chr);
	printf("******************NUM: %d\n",num);
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
void filterBlocks(){

// Duplicates in Genome A
	for(int i=0; i<BLOCK_NUM; i++){

		//std::cout<<"CHECK"<<std::endl;

	//	std::cout<<i<<"\n";
		int check =0;
		for(int j=i+1; j< BLOCK_NUM && check==0; j++){

	//		std::cout<<"CHECK1"<<std::endl;

	//		std::cout << blocks[i].achr<<"\n";

//			std::cout<<"CHECKsdvgasdvas1"<<std::endl;

			if(strcmp(blocks[i].achr, blocks[j].achr)==0){

				if(blocks[i].astart > blocks[j].astart){
					printf("Error in input file. Positions in the first genome need to be sorted and larger \n"
							"i.astart: %d \t j.astart: %d \n", blocks[i].astart, blocks[j].astart);
					exit(1);
				}
				if(blocks[i].astart == blocks[j].astart){
					if(blocks[i].aend > blocks[j].aend){
						blocks[j].state = DUP;
					}
					else{
						if(blocks[i].aend < blocks[j].aend){
						blocks[i].state = DUP;
						check = 1;
						}
						else{
							if(blocks[i].iden > blocks[j].iden){
								blocks[j].state = DUP;
							}
							else{
								blocks[i].state = DUP;
								check = 1;
							}
						}
					}
				}
				else{
					if(blocks[i].aend > blocks[j].astart){
						if(blocks[i].aend >= blocks[j].aend){
							blocks[j].state = DUP;
						}
						else{
							if(abs(blocks[i].astart - blocks[j].astart) < 50 && abs(blocks[i].aend - blocks[j].aend) < 50){
								if(blocks[i].state != DUP) blocks[i].state = TEMP_DUP;
								if(blocks[j].state != DUP) blocks[j].state = TEMP_DUP;
							}
							else{
								if(abs(blocks[i].astart - blocks[j].astart) > 50 && abs(blocks[i].aend - blocks[j].aend) < 50){
									blocks[j].state = DUP;
								}
								else{
									if(abs(blocks[i].astart - blocks[j].astart) < 50 && abs(blocks[i].aend - blocks[j].aend) > 50){
										blocks[i].state = DUP;
										check = 1;
									}
								}
							}
						}
					}
				}
			}
		}
	}

//	Duplicates in Genome B
	for(int i=0; i< BLOCK_NUM; i++){
	//	std::cout<<"CHECK"<<std::endl;
		if(blocks[i].state == DUP) continue;

		int check = 0;
	//	std::cout<<"CHECK"<<"\n";

		for(int j = i+1;j<BLOCK_NUM && check == 0; j++){
	//		std::cout<<"CHECK"<<"\n";

			if(blocks[j].state == DUP) continue;

			if(strcmp(blocks[i].bchr, blocks[j].bchr) == 0){
				if(blocks[i].bstart > blocks[j].bstart){
					if(blocks[i].bstart >= blocks[j].bend) continue;
					if (blocks[i].bend <= blocks[j].bend){
						blocks[i].state = DUP;
						check = 1;
					}
					else{
						if(abs(blocks[i].bstart - blocks[j].bstart) < 50 && abs(blocks[i].bend - blocks[j].bend) < 50){
							if(blocks[i].state != DUP) blocks[i].state = TEMP_DUP;
							if(blocks[j].state != DUP) blocks[j].state = TEMP_DUP;
						}
						else{
							if(abs(blocks[i].bstart - blocks[j].astart) > 50 && abs(blocks[i].aend - blocks[j].aend) < 50){
								blocks[i].state = DUP;
								check = 1;
							}
							else{
								if(abs(blocks[i].astart - blocks[j].astart) < 50 && abs(blocks[i].aend - blocks[j].aend) > 50){
									blocks[j].state = DUP;
								}
							}
						}
					}
				}
				else{
					if(blocks[i].bstart == blocks[j].bstart){
						if(blocks[i].bend > blocks[j].bend){
							blocks[j].state = DUP;
						}
						else{
							if(blocks[i].bend < blocks[j].bend){
								blocks[i].state = DUP;
								check = 1;
							}
							else{
								if(blocks[i].iden > blocks[j].iden){
									blocks[j].state = DUP;
								}
								else{
									blocks[i].state = DUP;
									check = 1;
								}
							}
						}
					}
					else{
						if(blocks[i].bend > blocks[j].bstart){
							if(blocks[i].bend >= blocks[j].bend){
								blocks[j].state = DUP;
							}
							else{
								if(abs(blocks[i].bstart - blocks[j].bstart) < 50 && abs(blocks[i].bend - blocks[j].bend) < 50){
									if(blocks[i].state != DUP) blocks[i].state = TEMP_DUP;
									if(blocks[j].state != DUP) blocks[j].state = TEMP_DUP;
								}
								else{
									if(abs(blocks[i].bstart - blocks[j].bstart) > 50 && abs(blocks[i].bend - blocks[j].bend) < 50){
										blocks[j].state = DUP;
									}
									else{
										if(abs(blocks[i].bstart - blocks[j].bstart) < 50 && abs(blocks[i].bend - blocks[j].bend) > 50){
											blocks[i].state = DUP;
											check = 1;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

}