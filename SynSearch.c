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
#include <assert.h>

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

    printf("parse ctx\n");

    parseCTX();
    printf("%s %d ***\n", blocks[58].achr, blocks[58].astart);

   // printf("lol\n");
//printf("%d\n", CHROMOSOME_NUM);

    for (int chr = 0; chr < CHROMOSOME_NUM; chr++) {

    //	printf("%d\n",chr);
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

    //    printf("%s %d 2\n", blocks[58].achr, blocks[58].astart);

   		for(int i=0; i < BLOCK_NUM; i++){
    		if(strcmp(blocks[i].achr, CHROMOSOME[chr])== 0){
    			chromo[chromo_index] = blocks[i];
    			strcpy(chromo[chromo_index].achr, blocks[i].achr);
    			strcpy(chromo[chromo_index].bchr, blocks[i].bchr);

    			chromo_index++;
    		}
    	}

   	 //   printf("%s %d 3\n", chromo[59].achr, chromo[59].astart);


   		SYNPATH synPath = parseSYN(chromo, CHROMOSOME[chr], num);

   	//    printf("%s %d 4\n", chromo[59].achr, chromo[59].astart);


   		printSynPath(chromo, synPath);

   	   // printf("%s %d 5\n", chromo[59].achr, chromo[59].astart);


   		parseITX(chromo, CHROMOSOME[chr], num);

   		parseINV(chromo, CHROMOSOME[chr], num);

   		groupITX(chromo, CHROMOSOME[chr], num);
   		int a =15;
   		assert(a >= 10);

   		free(chromo);
    	chromo = NULL;
    	//free(maxWeightPath);
    	//free(longestInverted);
    	free(synPath.maxWeightPath);
    }

    printf("\n\n FINISHED SUCCESSFULLY");

    fclose(itxOutFile);
	return EXIT_SUCCESS;
}

