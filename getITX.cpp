/*
 * getITX.c
 *
 *  Created on: Feb 2, 2017
 *      Author: goel
 */

#include "getITX.h"


void parseITX(std::vector<BLOCK> &chromo, char chr[], int num) {
	int i;
	for (i = 1; i < num-1; i++) {
		if (chromo[i].bchr.compare(chr) == 0 && chromo[i].state != CTX && chromo[i].state != SYN) { //&& blocks[i].state != SYN_IN_INV && blocks[i].state != INV) {
			if(chromo[i].dir == 1){
				chromo[i].state = ITX;
			}
			else if (chromo[i].dir == -1){
				chromo[i].state = INV_ITX;
			}
		}
	}
}



void groupITX(std::vector<BLOCK> &chromo, char chr[], int num){
	int i, in = -1;
	int state;


	for (i = 1; i < num-1; i++) {

		if (chromo[i].bchr.compare(chr) == 0 && (chromo[i].state == ITX || chromo[i].state == INV_ITX || chromo[i].state == ITX_IN_INV)) {

			if (in == -1) {
				in = i;
				state = chromo[i].state;
			}
			else {
				// There is a itx already, do they belong together?
				// If yes, do nothing
				// If no, print old set
				int together = 1;
				if (chromo[i].bchr.compare(chromo[i-1].bchr) != 0 || chromo[i].dir != chromo[i-1].dir) {
					// not on same chromosome or different direction?
					together = 0;

				}
				else { // on same chr and same direction, but next to each other also on B genome?
					if ((chromo[i].dir == 1 && chromo[i].leftBNeighbor != i-1) ||
							(chromo[i].dir == -1 && chromo[i].rightBNeighbor != i-1)) {
						together = 0;
					}
					else{
						if(chromo[i].state != state){
							together = 0;
						}
					}

				}
				if (together == 0) {
					writeITX(chromo, in, i-1);
					in = i;
				}
			}
		}
		else {
			// Write if a ITX was just passed
			if (in != -1) {
				writeITX(chromo, in, i-1);
			}
			in = -1;
		}
	}

	if (in != -1) {
		writeITX(chromo,in, i-1);
	}
}


void writeITX(std::vector<BLOCK> &chromo, int a, int b) {

	//Header
	switch(chromo[a].state){
	case 5://ITX
		fprintf(itxOutFile,"#ITX\t");
		//printf("chr : %s \t chr: %d\n", chromo[a].achr, chromo[a].achr);
		fprintf(itxOutFile, "%s %d %d - ", chromo[a].achr.c_str(), chromo[a].astart, chromo[b].aend);
		fprintf(itxOutFile, "%s %d %d\n", chromo[a].bchr.c_str(), chromo[b].bstart, chromo[b].bend);
		break;

	case 6: //Inverted ITX
		fprintf(itxOutFile,"#INV_ITX\t");
		fprintf(itxOutFile, "%s %d %d - ", chromo[a].achr.c_str(), chromo[a].astart, chromo[b].aend);
		fprintf(itxOutFile, "%s %d %d\n", chromo[a].bchr.c_str(), chromo[b].bstart, chromo[a].bend);
		break;

	case 7: //Inverted ITX
		fprintf(itxOutFile,"#ITX_IN_INV\t");
		fprintf(itxOutFile, "%s %d %d - ", chromo[a].achr.c_str(), chromo[a].astart, chromo[b].aend);
		fprintf(itxOutFile, "%s %d %d\n", chromo[a].bchr.c_str(), chromo[b].bstart, chromo[a].bend);
		break;
	}

for (int i = a; i <= b; i++) {
		writeBlock(itxOutFile, chromo[i]);
	}

}
