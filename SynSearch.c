/*
 ============================================================================
 Name        : SynSearch.c
 Author      : Korbinian Schneeberger
 Version     : 0.1
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "SynSearch.h"

/*
 * Missing:
 * Handle inversions at end of chr properly (perhaps introduce syntenic pseudo blocks at each chr end)
 * CTX parsing could already profit from B genome edges
 * Annotate Syn_in_inversions not in the syntenic path (?)
 * Allow for generic output file names
 *
 * Define gaps between syntenics blocks: different types of complex regions
 * It is not clear how MUMmer handles duplications: simulate and test. Delta-filter with parameter -1 filters for one-to-one alignments (not to confuse with -l)
 * Nucmer's --simplify might remove many duplications, which is necessary for Synsearch (so far)
 */

int main(int argc, char *argv[]) {
       init(argc, argv);
    printf("init\n");

    //fills some global variables, most importantly it fills blocks
    readInputFile();
    printf("read file\n");

    // First: parse out cross chromosomal translocations
    parseCTX();
    printf("parse ctx\n");

    for (int chr = 0; chr < CHROMOSOME_NUM; chr++) {

    	//Connect B genome
		setEdgesBGenome(CHROMOSOME[chr]);

		/////////////////////////////////
		// Find syntenic path

		// build all possible paths
		setEdges(CHROMOSOME[chr]);

		// calculate cummulative weigths for each path
		// while only the heaviest path is followed
		setPathWeights(CHROMOSOME[chr]);

		// fills maxWeightPath[] with the co-linear blocks
		// and sets their status to SYN
		backtraceSynPath();

		// print path block and combine runs of blocks
		printSynPath();

		maxWeight = 0;

		printf("parse syn regions\n");

		// Find inversions
		// This is not optimally solved. Inv at the end of chr
		// are not found, and ITX within the inversion are also
		// not reported. (In fact this would be higher level #3).
              printf("%s\n",CHROMOSOME[chr]);

 		writeInversions(CHROMOSOME[chr]);
		printf("write inversions\n");

	    // write out ITX (all the rest, but need to be combined similar to the CTX)
	    // does includes ITX(inverted)
	    parseITX(CHROMOSOME[chr]);

    }

       printf("\n\n FINISHED SUCCESSFULLY");
	return EXIT_SUCCESS;
}

void parseITX(char chr[]) {
	int i, in = -1;

	for (i = 0; i < BLOCK_NUM; i++) {

		if (strcmp(blocks[i].bchr, chr) == 0 && blocks[i].state != CTX && blocks[i].state != SYN && blocks[i].state != SYN_IN_INV && blocks[i].state != INV) {
			blocks[i].state = ITX;

			//if (blocks[i].astart == 11394535 || blocks[i].astart == 11408511) {
			//	printBlock(i);
			//}

			if (in == -1) {
				in = i;
			}
			else {
				// There is a ctx already, do they belong together?
				// If yes, do nothing
				// If no, print old set
				int together = 1;
				if (strcmp(blocks[i].bchr, blocks[i-1].bchr) != 0 || blocks[i].dir != blocks[i-1].dir) {
					// not on same chromosome or different direction?
					together = 0;

				}
				else { // on same chr and same direction, but next to each other also on B genome?
					if ((blocks[i].dir == 1 && blocks[i].leftBNeighbor != i-1) ||
							(blocks[i].dir == -1 && blocks[i].rightBNeighbor != i-1)) {
						together = 0;
					}
				}

				if (together == 0) {
					writeITX(in, i-1);
					in = i;
				}
			}

		}
		else {
			// Write if a CTX was just passed
			if (in != -1) {
				writeITX(in, i-1);
			}
			in = -1;
		}
	}

	if (in != -1) {
		writeCTX(in, i-2);
	}

}

void writeITX(int a, int b) {

	//Header
	fprintf(itxOutFile, "#ITX");
	if (blocks[a].dir == -1) {
		fprintf(itxOutFile, "(inverted) ");
	}
	else {
		fprintf(itxOutFile, " ");
	}
	fprintf(itxOutFile, "%s %d %d - ", blocks[a].achr, blocks[a].astart, blocks[b].aend);
	if (blocks[a].dir == 1) {
		fprintf(itxOutFile, "%s %d %d\n", blocks[a].bchr, blocks[a].bstart, blocks[b].bend);
	}
	else {
		fprintf(itxOutFile, "%s %d %d\n", blocks[a].bchr, blocks[b].bstart, blocks[a].bend);
	}

	//all blocks
	for (int i = a; i <= b; i++) {
		writeBlock(itxOutFile, blocks[i]);
	}

}


void writeInversions(char chr[]) {

	for (int i = 0; i < BLOCK_NUM; i++) {
		//printf("%d (%d)\n", i, BLOCK_NUM);
		if (strcmp(blocks[i].bchr, chr) == 0 && blocks[i].state != CTX
			//&& blocks[i].astart != 16047157

		    ){
			if (blocks[i].dir == -1) {
				// now decide whether the inversion is translocated or not
				// for this first get first
				//if (blocks[i].rightNeighbor )

				/////////////////////////////////////////////////////
				// Set right B neighbor
				/////////////////////////////////////////////////////
				int rightB = blocks[i].rightBNeighbor;
				while (rightB != -1 && blocks[rightB].state != SYN) {
					rightB = blocks[rightB].rightBNeighbor;
				}
				//printf("rightB %d\n", rightB);

				/////////////////////////////////////////////////////
				// Set left A neighbor
				/////////////////////////////////////////////////////
				int leftA = i - 1;
				while (leftA != -1 && blocks[leftA].state != SYN) {
					leftA--;
				}
				//printf("leftA %d\n", leftA);

				if (!(leftA == -1 && rightB == -1)) { //Special case: whole chr inverted, means no inversion as then the alignment dir wrong

					if (rightB == -1 || leftA == -1) {
				//		printf("WARNING: Inv at the end of chr: not solved. Will be reported as ITX.\n");  NEED TO UNCOMMENT
					}
					else {

						if (rightB > i) { // otherwise running back is senseless and no inversion is found

							//Run back on A genome until an inversion is found
							int invEndA = rightB;
							/** if (rightB == -1) { // inv at the end of the chromosome
								invEndA = BLOCK_NUM-1;
							}*/
							while (blocks[invEndA].dir == 1 || blocks[invEndA].state == CTX || blocks[invEndA].state == SYN) {
								invEndA--;
							}

							//printf("invEndA %d\n", invEndA);

							//Run left until next syntenic block
							int leftB = blocks[invEndA].leftBNeighbor;
							/** MISSING SOLUTION END OF CHR */
							while (leftB != -1 && blocks[leftB].state != SYN) {
								leftB = blocks[leftB].leftBNeighbor;
							}
							//printf("leftB %d\n", leftB);


							if (leftB == -1) {
								printf("WARNING: Inv at the end of chr: not solved. Will be reported as ITX.\n");
							}
							else {

								// Run back on A genome until an inversion is found
								int invBeginA = leftB;
								/**if (leftB == -1) {
									invBeginA = 0;
								}*/

						//		printf("random check\n");

								while (blocks[invBeginA].dir != -1 || blocks[invBeginA].state == CTX || blocks[invBeginA].state == SYN) {
									invBeginA++;
									if(invBeginA > BLOCK_NUM){
										break;
									}
								}
						//		printf("random check2\n");

								//printf("invBeginA %d\n", invBeginA);

								// If the end of the circuit equals the beginning: Inversion found
								// This can include rearranged inversions in the inversion.
								// Actually the best would be to take inversions and run SnyPathFinder within them (respecting the inversion)
								// This would give the optimal inversion and report the ITX in the inversion
								if (invBeginA == i) {
									//printf("INVERSION! invBeginBlock %d invEndBlock %d\n", invBeginA, invEndA);
									fprintf(invOutFile, "#INV ");
									fprintf(invOutFile, "%s %d %d - ", blocks[invBeginA].achr, blocks[invBeginA].astart, blocks[invEndA].aend);
									fprintf(invOutFile, "%s %d %d\n", blocks[invEndA].bchr, blocks[invEndA].bstart, blocks[invBeginA].bend);

									for (int j = invBeginA; j <= invEndA; j++) {
										if (blocks[j].state == SYN) {
											blocks[j].state = SYN_IN_INV;
										}
										else {
											if (blocks[j].state == CTX) {
												//do nothing
											}
											else {
												if (blocks[j].dir == -1) {
													blocks[j].state = INV;
													writeBlock(invOutFile, blocks[j]);
												}
											}
										}
									}
									i = invEndA;
								}

							}
						}
					}
				}

			}
		}
	}
       //printf("****Finished Writing Inversions****\n");

}


void setEdgesBGenome(char chr[]) {

	int first = -1, last = -1;


	for (int i = 0; i < BLOCK_NUM; i++) {
		blocks[i].leftBNeighbor = -1;
		blocks[i].rightBNeighbor = -1;
		if (strcmp(blocks[i].bchr, chr) == 0 && blocks[i].state != CTX){ // only one chr at a time and only on one chr
			// This is not a good sorting algorithm,
			// but as the list is almost sorted
			// this should work.

			//initiate list:
			if (last == -1) {
				first = i;
				last = i;
			}
 			// if there is at least one element:
			else {
				int current = last;
				//Move from end into list to find correct position:
				//printf("i %d   last %d current %d current.leftN %d  current.leftN[5489] %d  current.leftN[858] %d\n", i, last, current, blocks[current].leftBNeighbor, blocks[5489].leftBNeighbor, blocks[858].leftBNeighbor);
				while (current != -1 && blocks[current].bstart > blocks[i].bstart) {
					current = blocks[current].leftBNeighbor;
					//printf("c: %d\n", current);
				}
				if (current == -1) { // i sits at the beginning
					blocks[i].rightBNeighbor = first;
					blocks[first].leftBNeighbor = i;
					first = i;
				}
				else {
					blocks[i].leftBNeighbor = current;
					blocks[i].rightBNeighbor = blocks[current].rightBNeighbor;

					blocks[blocks[i].leftBNeighbor].rightBNeighbor = i;
					if (blocks[i].rightBNeighbor != -1) {
						blocks[blocks[i].rightBNeighbor].leftBNeighbor = i;
					}
					else {
						last = i;
					}
				}
			}
		}
	}
}

void parseCTX() {

	int i, in = -1;

	//Parse through all blocks to find translocations between chromosomes
	//those blocks that come from the same translocation, should be reported together
	for (i = 0; i < BLOCK_NUM; i++) {
		if (strcmp(blocks[i].achr, blocks[i].bchr) != 0) { // yes the two blocks are assigned to different chromosomes
			blocks[i].state = CTX;

			if (in == -1) {
				in = i;
			}
			else {
				// There is a ctx already, do they belong together?
				// If yes, do nothing
				// If no, print old set
				int together = 1;
				if (strcmp(blocks[i].bchr, blocks[i-1].bchr) != 0 || blocks[i].dir != blocks[i-1].dir) {
					// not on same chr or diff direction?
					together = 0;
				}
				else { // on same chr, but next to each other?
					for (int j = 0; j < BLOCK_NUM; j++) {
						if (j != i && j != i-1) {
							// Is block j between i and i-1 on the b-genome?
							if ((blocks[j].bstart > blocks[i].bstart && blocks[j].bstart < blocks[i-1].bstart) // rev
								|| (blocks[j].bstart < blocks[i].bstart && blocks[j].bstart > blocks[i-1].bstart)) {
									together = 0;
							}
						}
					}
				}

				if (together == 0) {
					writeCTX(in, i-1);
					in = i;
				}
			}

		}
		else {
			// Write if a CTX was just passed
			if (in != -1) {
				writeCTX(in, i-1);
			}
			in = -1;
		}
	}

	if (in != -1) {
		writeCTX(in, i-2);
	}

}

void writeCTX(int a, int b) {

	//Header
	fprintf(ctxOutFile, "#CTX");
	if (blocks[a].dir == -1) {
		fprintf(ctxOutFile, "(inverted) ");
	}
	else {
		fprintf(ctxOutFile, " ");
	}
	fprintf(ctxOutFile, "%s %d %d - ", blocks[a].achr, blocks[a].astart, blocks[b].aend);
	if (blocks[a].dir == 1) {
		fprintf(ctxOutFile, "%s %d %d\n", blocks[a].bchr, blocks[a].bstart, blocks[b].bend);
	}
	else {
		fprintf(ctxOutFile, "%s %d %d\n", blocks[a].bchr, blocks[b].bstart, blocks[a].bend);
	}

	//all blocks
	for (int i = a; i <= b; i++) {
		writeBlock(ctxOutFile, blocks[i]);
	}

}

void backtraceSynPath() {

	int maxWeightPathReverse[MAX_BLOCK_NUM];
	int length = 0;

	int cblock = maxWeightBlock;

	while (blocks[cblock].inEdgeNum > 0) {
		maxWeightPathReverse[length] = cblock;
		cblock = blocks[cblock].maxInEdge;
		length++;
	}

	maxWeightPathReverse[length] = cblock;
	length++;

	for (int i = length-1; i >= 0; i--) {
		maxWeightPath[length-i-1] = maxWeightPathReverse[i];
		blocks[maxWeightPath[length-i-1]].state = SYN;
	}

	maxWeightPathLength = length;
}

void printSynPath() {

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
			if (i-1 != maxWeightPath[p-1] || blocks[i].leftBNeighbor != maxWeightPath[p-1]) { // new region
				fprintf(synOutFile, "#SYN ");
				fprintf(synOutFile, "%s %d %d - ", blocks[maxWeightPath[s]].achr, blocks[maxWeightPath[s]].astart, blocks[maxWeightPath[p-1]].aend);
				fprintf(synOutFile, "%s %d %d\n", blocks[maxWeightPath[s]].bchr, blocks[maxWeightPath[s]].bstart, blocks[maxWeightPath[p-1]].bend);
				for (int j = s; j < p; j++) {
					writeBlock(synOutFile, blocks[maxWeightPath[j]]);
				}
				s = p;
				in = 1;
			}
		}
	}

	if (in == 1) {
		int e = maxWeightPathLength-1;
		fprintf(synOutFile, "#SYN ");
		fprintf(synOutFile, "%s %d %d - ", blocks[maxWeightPath[s]].achr, blocks[maxWeightPath[s]].astart, blocks[maxWeightPath[e]].aend);
		fprintf(synOutFile, "%s %d %d\n", blocks[maxWeightPath[s]].bchr, blocks[maxWeightPath[s]].bstart, blocks[maxWeightPath[e]].bend);
		for (int j = s; j < maxWeightPathLength; j++) {
			writeBlock(synOutFile, blocks[maxWeightPath[j]]);
		}
	}
}

void setPathWeights(char chr[]) {

	for (int i = 0; i < BLOCK_NUM; i++) {
		if (strcmp(blocks[i].achr, chr) == 0 && blocks[i].dir == 1 && blocks[i].state != CTX){ // only one chromosome at a time and only fwd blocks
			blocks[i].weight = blocks[i].alen;
			for (int j = 0; j < blocks[i].inEdgeNum; j++) {
				if (blocks[i].weight < blocks[blocks[i].inEdge[j]].weight + blocks[i].alen) {
					blocks[i].weight = blocks[blocks[i].inEdge[j]].weight + blocks[i].alen;
					blocks[i].maxInEdge = blocks[i].inEdge[j];
					if (blocks[i].weight > maxWeight) {
						maxWeightBlock = i;
						maxWeight = blocks[i].weight;
					}
				}
			}
		}
	}

	//printf("last block of best path %d\n", maxWeightBlock);
	//printf("max weight %d\n", maxWeight);
}

void setEdges(char chr[]) {

	for (int i = BLOCK_NUM-1; i >= 0; i--) {
		if (strcmp(blocks[i].achr, chr) == 0 && blocks[i].dir == 1 && blocks[i].state != CTX){ // only one chromosome at a time
			int setEdge = 0;
			int j = i+1;
			blocks[i].flagBadOut = 1;
			while (setEdge == 0 && j < BLOCK_NUM) { // block j needs to exist
				if (strcmp(blocks[j].achr, chr) == 0 && blocks[j].dir == 1) {
					if (testSynteny(i, j)) {
						blocks[i].outEdge[blocks[i].outEdgeNum] = j;
						blocks[i].outEdgeNum++;
						blocks[j].inEdge[blocks[j].inEdgeNum] = i;
						blocks[j].inEdgeNum++;
						// set badOut to 0 if it is the next block that is linked
						if (i+1 == j) {
							blocks[i].flagBadOut = 0;
						}
						// Only if the linked block is not badOut stop
						// This is relevant if multiple translocations follow each other: multiple outEdges will be set for the block
						if (blocks[j].flagBadOut == 0) {
							setEdge = 1;
						}
					}
				}
				j++;
			}
		}
	}

}

int testSynteny(int i, int j) {
	if (blocks[i].astart < blocks[j].astart && blocks[i].bstart < blocks[j].bstart && strcmp(blocks[i].achr, blocks[j].bchr)==0) {
		return 1;
	} else {
		return 0;
	}
}

void readInputFile() {

  if ((inputFile = fopen(inputFileName, "r")) == NULL) {
	  printf("Cannot open input file\n");
    exit(1);
  }

  int i1, i2, i3, i4, i5, i6, i8, i9;
  float f7, f10, f11;
  char c12[4096], c13[4096];
  int ret=0, adir, bdir;
  int c = 0;
  int line = 0;
  char currchr[4096];
  int currpos = 0;

  //read in file

  while((ret = fscanf(inputFile, "%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%s\t%s\n",
		  &i1, &i2, &i3, &i4, &i5, &i6, &f7, &i8, &i9, &c12[0], &c13[0])) == 11) {

	  line++;

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

	  if (adir == bdir) {
		  blocks[c].dir = 1;
	  } else {
		  blocks[c].dir = -1;
	  }

	  //record chromosome/scaffold/contig of that block in both genomes
	  strcpy(blocks[c].achr, c12);
	  strcpy(blocks[c].bchr, c13);

	  // increment block counter
	  c++;

	  // collect chromosome identifier:
	  if (strcmp(c12, currchr) == 0) { // same chr
		  if (i1 < currpos) {
			  printf("Error in input file. Positions in the first genome need to be sorted and larger 0 (line %d)\n", line);
			  exit(1);
		  }
		  currpos = i1;
	  }
	  else { // different chromosome
		  for (int i = 0; i < CHROMOSOME_NUM; i++) { // there already?
				 if (strcmp(CHROMOSOME[i], c12) == 0) {
					 printf("Error in input file. Chromosomes are not sorted (line %d).\n", line);      // should add an exit command?
				 }
		  }
		  strcpy(CHROMOSOME[CHROMOSOME_NUM], c12);
		  CHROMOSOME_NUM++;
		  strcpy(currchr, c12);
		  currpos = i1;
	  }
  }

  BLOCK_NUM = c;

  if (ret != EOF) {
      printf("Input file error\n");
      exit(1);
  }

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

void printBlocks() {
	for (int i = 0; i < BLOCK_NUM; i++) {
		printBlock(i);
	}
}

void printBlock(int i) {
	printf("################################################\n");
	printf("Block %d\n", i);
	printf(" Genome A: %s   %d   ---   %d\n", blocks[i].achr, blocks[i].astart, blocks[i].aend);
	printf(" Genome B: %s   %d   ---   %d\n", blocks[i].bchr, blocks[i].bstart, blocks[i].bend);
	printf(" Direction: %d\n", blocks[i].dir);
	printf(" State: %d\n", blocks[i].state);
	printf("\n");
	printf("Badout?: %d\n", blocks[i].flagBadOut);
	printf("Out-Edge num: %d\n", blocks[i].outEdgeNum);
	printf("Out blocks:");
	for (int j = 0; j < blocks[i].outEdgeNum; j++) { printf(" %d", blocks[i].outEdge[j]); }
	printf("\n");
	printf("\n");
	printf("In-Edge num: %d\n", blocks[i].inEdgeNum);
	printf("In blocks:");
	for (int j = 0; j < blocks[i].inEdgeNum; j++) { printf(" %d", blocks[i].inEdge[j]); }
	printf("\n");
	printf("\n");
	printf("Weight: %d\n", blocks[i].weight);
	printf("Max in Edge: %d\n", blocks[i].maxInEdge);
	printf("\n");
	printf("Left neighbor: %d\n", blocks[i].leftBNeighbor);
	printf("Right neighbor: %d\n", blocks[i].rightBNeighbor);
}

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
  printf("SynSearch, version 0.02, 09.02.2017\n");
}

