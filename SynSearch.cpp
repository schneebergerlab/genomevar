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
#include<iostream>
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

using namespace std;

int main(int argc, char *argv[]) {


	init(argc, argv);

    //fills some global variables, most importantly it fills blocks
  //  readInputFile();
    printf("read file\n");

    parseDUP(blocks, mblocks);


    int countUNI=0, countDUP=0, countRED = 0;

    for(int i=0; i<BLOCK_NUM; ++i){
    	if(blocks[i].state == UNI) ++countUNI;
    	if(blocks[i].state == DUP) ++countDUP;
    	if(blocks[i].state == RED) ++countRED;
    }

    std::cout<<countUNI<<"\t"<<countDUP<<"\t"<<countRED<<"\n";

    countUNI=0;
    countDUP=0;

    for(int i=0; i<mBLOCK_NUM; ++i){
       	if(mblocks[i].state == UNI) ++countUNI;
       	if(mblocks[i].state == DUP) ++countDUP;
       }

    std::cout<<countUNI<<"\t"<<countDUP<<"\n";

    // First: parse out cross chromosomal translocations

    printf("parse ctx\n");
    parseCTX();

int dup_count =0;


//=0;
for(unsigned i = 0; i < blocks.size();++i){
 	if(blocks[i].state == DUP){
		dup_count++;
	}

}
   for (int chr = 0; chr < CHROMOSOME_NUM; chr++) {
    	std::vector<int> indices;
    	for ( int i =0; i < BLOCK_NUM; i++){
    		if(strcmp(blocks[i].achr, CHROMOSOME[chr])== 0){
    			indices.push_back(i);
    		}
    	}


    	std::vector<BLOCK> chromo(indices.size()); // = (BLOCK *) calloc(num, sizeof(BLOCK));

    //	int chromo_index=1;

   		for(std::vector<int>::iterator it = indices.begin();it!=indices.end();++it){
   			chromo[it-indices.begin()] = blocks[*it];
   		}

   		BLOCK new_block;

   		new_block.state = STER;
   		chromo.insert(chromo.begin(), new_block);

   		new_block.state = ETER;
   		chromo.push_back(new_block);

   		SYNPATH synPath = parseSYN(chromo, CHROMOSOME[chr], chromo.size());

   		printSynPath(chromo, synPath);


   		parseITX(chromo, CHROMOSOME[chr],  chromo.size());

   		parseINV(chromo, CHROMOSOME[chr],  chromo.size());

   		groupITX(chromo, CHROMOSOME[chr],  chromo.size());

   		//free(chromo);
 //   	chromo = NULL;
    	//free(maxWeightPath);
    	//free(longestInverted);
    	free(synPath.maxWeightPath);


    }

    printf("\n\n FINISHED SUCCESSFULLY");

    fclose(itxOutFile);
    fclose(synOutFile);
    fclose(invOutFile);
    fclose(ctxOutFile);
    fclose(dupOutFile);
	return EXIT_SUCCESS;
	/*free(mBlocks);
	return(0);*/

}

