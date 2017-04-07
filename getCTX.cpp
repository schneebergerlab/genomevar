/*
 * getCTX.c
 *
 *  Created on: Feb 2, 2017
 *      Author: goel
 */
#include "init.h"
#include "getCTX.h"

void parseCTX(std::vector<BLOCK> &blocks) {


    int bSize = blocks.size();

    std::cout<<"bsize: "<<bSize<<"\n";
	int i, in = -1;
	//Parse through all blocks to find trans-locations between chromosomes
	//those blocks that come from the same translocation, should be reported together
	for (i = 0; i < bSize; i++) {
		if (blocks[i].achr.compare(blocks[i].bchr) != 0) { // yes the two blocks are assigned to different chromosomes
            //std::cout<<"LOL1\n";
            blocks[i].state = CTX;
			if (in == -1) {
				in = i;
			}
			else{
				// There is a ctx already, do they belong together?
				// If yes, do nothing
				// If no, print old set
				int together = 1;
				if (blocks[i].bchr.compare(blocks[i-1].bchr) != 0 || blocks[i].dir != blocks[i-1].dir) {
					// not on same chr or diff direction?
					together = 0;
				}
				else { // on same chr, but next to each other?
					for (int j = 0; j < bSize; j++) {
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
					writeCTX(blocks, in, i-1);
					in = i;
				}
			}

		}
		else {
			// Write if a CTX was just passed
			if (in != -1) {
				writeCTX(blocks, in, i-1);
			}
			in = -1;
		}
	}
	if (in != -1) {
		writeCTX(blocks, in, i-1);
	}

	std::cout<<"DONE\n";
}

void writeCTX(std::vector<BLOCK> &blocks, int a, int b) {
	//Header
	fprintf(ctxOutFile, "#CTX");
	if (blocks[a].dir == -1) {
		fprintf(ctxOutFile, "(inverted) ");
	}
	else {
		fprintf(ctxOutFile, " ");
	}

	fprintf(ctxOutFile, "%s %d %d - ", blocks[a].achr.c_str(), blocks[a].astart, blocks[b].aend);

	if (blocks[a].dir == 1) {
		fprintf(ctxOutFile, "%s %d %d\n", blocks[a].bchr.c_str(), blocks[a].bstart, blocks[b].bend);
	}
	else {
		fprintf(ctxOutFile, "%s %d %d\n", blocks[a].bchr.c_str(), blocks[b].bstart, blocks[a].bend);
	}

	//all blocks
	for (int i = a; i <= b; i++) {
		writeBlock(ctxOutFile, blocks[i]);
	}
}
