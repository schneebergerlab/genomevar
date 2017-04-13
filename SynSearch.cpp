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
    std::vector<BLOCK> uniBlocks;
    std::vector<BLOCK> dupBlocks;

	init(argc, argv);
	//estimateThreshold(blocks);

    printf("read file\n");

    if(argc == 5 and strcmp(argv[1], "-1")==0 and strcmp(argv[3],"-m")==0){
        FILTEREDDATA fData;
        parseDUP(blocks, mblocks, 50, fData);
        uniBlocks = fData.uniBlocks;
        dupBlocks = fData.dupBlocks;
	}
	else{
        uniBlocks = blocks;
	}


    int uniSize = uniBlocks.size();
    parseCTX(uniBlocks);
    for (int chr = 0; chr < CHROMOSOME_NUM; chr++) {
    	std::vector<int> indices;
    	for ( int i =0; i < uniSize; ++i){
    		if(uniBlocks[i].achr.compare(CHROMOSOME[chr])== 0){
    			indices.push_back(i);
    		}
    	}

    	std::vector<BLOCK> chromo(indices.size());
   		for(std::vector<int>::iterator it = indices.begin();it!=indices.end();++it){
   			chromo[it-indices.begin()] = uniBlocks[*it];
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
    dupOutFile.close();
	return EXIT_SUCCESS;
	/*free(mBlocks);
	return(0);*/

}

