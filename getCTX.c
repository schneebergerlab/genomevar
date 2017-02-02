/*
 * getCTX.c
 *
 *  Created on: Feb 2, 2017
 *      Author: goel
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "init.h"
#include "getCTX.h"

void parseCTX() {

	int i, in = -1;
	//Parse through all blocks to find trans-locations between chromosomes
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

	printf("FINISHED CTX\n");

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


