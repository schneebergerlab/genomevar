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
 //   readInputFile();
    printf("read file\n");

    // First: parse out cross chromosomal translocations
    parseCTX();

    printf("parse ctx\n");

    for (int chr = 0; chr < CHROMOSOME_NUM; chr++) {
    	int num=0;
    	for ( int i =0; i < BLOCK_NUM; i++){
    		if(strcmp(blocks[i].achr, CHROMOSOME[chr])== 0){
    			num++;
    		}
    	}

    	num = num+2;				//Size of chromosone with add start and terminal blocks

    	BLOCK *chromo = (BLOCK *) calloc(num, sizeof(BLOCK));

    	if (chromo == 0)
    	{
    		printf("ERROR: Out of memory\n");
    		return 1;
    	}

    	int chromo_index=1;
    	chromo[0].state =  STER;
    	chromo[num-1].state = ETER;


   		for(int i=0; i < BLOCK_NUM; i++){
    		if(strcmp(blocks[i].achr, CHROMOSOME[chr])== 0){
    			chromo[chromo_index] = blocks[i];
    			chromo_index++;
    		}
    	}
   	    parseSYN(chromo, CHROMOSOME[chr], num);
/*
   		//printf("setting B edges\n");
   		setEdgesBGenome(CHROMOSOME[chr], num);

   		// build all possible paths
		setEdges(CHROMOSOME[chr], num);

		// calculate cummulative weigths for each path
		// while only the heaviest path is followed
		setPathWeights(CHROMOSOME[chr], num);

		backtraceSynPath(CHROMOSOME[chr], num);

		printSynPath(CHROMOSOME[chr], num);
		maxWeight = 0;

		parseITX(CHROMOSOME[chr], num);

 		writeInversions(CHROMOSOME[chr], num);


    	free(chromo);
    	free(maxWeightPath);
    	//free(longestInverted);*/
    }


       printf("\n\n FINISHED SUCCESSFULLY");
	return EXIT_SUCCESS;
}

