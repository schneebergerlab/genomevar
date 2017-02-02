/*
 * getITX.c
 *
 *  Created on: Feb 2, 2017
 *      Author: goel
 */

#include "getITX.h"

/*
void parseITX(char chr[], int num) {
	int i;
	for (i = 1; i < num-1; i++) {
		if (strcmp(chromo[i].bchr, chr) == 0 && chromo[i].state != CTX && chromo[i].state != SYN) { //&& blocks[i].state != SYN_IN_INV && blocks[i].state != INV) {
			if(chromo[i].dir == 1){
				chromo[i].state = ITX;
			}
			else if (chromo[i].dir == -1){
				chromo[i].state = ITX_IN;
			}
		}
	}
}
*/


/*
void groupITX(char chr[]){

//			//if (blocks[i].astart == 11394535 || blocks[i].astart == 11408511) {
//			//	printBlock(i);
//			//}
//
//			if (in == -1) {
//				in = i;
//			}
//			else {
//				// There is a ctx already, do they belong together?
//				// If yes, do nothing
//				// If no, print old set
//				int together = 1;
//				if (strcmp(blocks[i].bchr, blocks[i-1].bchr) != 0 || blocks[i].dir != blocks[i-1].dir) {
//					// not on same chromosome or different direction?
//					together = 0;
//
//				}
//				else { // on same chr and same direction, but next to each other also on B genome?
//					if ((blocks[i].dir == 1 && blocks[i].leftBNeighbor != i-1) ||
//							(blocks[i].dir == -1 && blocks[i].rightBNeighbor != i-1)) {
//						together = 0;
//					}
//				}
//
//				if (together == 0) {
//					writeITX(in, i-1);
//					in = i;
//				}
//			}
//
//		}
//		else {
//			// Write if a CTX was just passed
//			if (in != -1) {
//				writeITX(in, i-1);
//			}
//			in = -1;
//		}
//	}
//
//	if (in != -1) {
//		writeCTX(in, i-2);
//	}
//
}
*/

