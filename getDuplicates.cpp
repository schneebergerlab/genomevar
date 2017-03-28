/*
 * getDuplicates.cpp
 *
 *  Created on: Mar 16, 2017
 *      Author: goel
 */
#include "getDuplicates.h"
//#include <iostream>

FILTEREDDATA parseDUP(std::vector<BLOCK> &blocks,std::vector<BLOCK> &mBlocks, int threshold){

	std::cout<< "FINDING DUPLICATES\n";

	int size = blocks.size();
    int msize = mBlocks.size();


	filterBlocks(blocks, threshold);
	filterBlocks(mBlocks, blocks, threshold);

	std::vector<BLOCK> mUni;
	for(int i =0; i < msize;++i){
		if(mBlocks[i].state == UNI){
			mUni.push_back(mBlocks[i]);
		}
	}

	int mUniSize = mUni.size();
	filterBlocks(mUni, threshold);

//	int count = 0;
//	int countDUP = 0;
//	for(int i = 0; i < mUni.size(); ++i){
//		/*if(mUni[i].state == UNI){
//			count++;
//			continue;
//		}*/
//		if(mUni[i].state == DUP){
//		mUni.erase(mUni.begin()+i);
//		continue;}
//	}

	FILTEREDDATA fData;

	for(int i = 0; i < size; ++i){
		if(blocks[i].state == UNI){
			fData.uniBlocks.push_back(blocks[i]);
			}
			else if(blocks[i].state == DUP){
                fData.dupBlocks.push_back(blocks[i]);
			}
        }

    for(int i = 0; i < msize; ++i){
        if(mBlocks[i].state == DUP) fData.dupBlocks.push_back(mBlocks[i]);
    }



    for(int i = 0; i < mUniSize; ++i){
        if(mUni[i].state == UNI){
            fData.uniBlocks.push_back(mUni[i]);
        }
        else if (mUni[i].state == DUP){
            fData.dupBlocks.push_back(mUni[i]);
        }
    }



    std::stable_sort(fData.uniBlocks.begin(), fData.uniBlocks.end(), []( const BLOCK &lhs, const BLOCK &rhs){
   // std::cout<<"lhs.start: "<<lhs.astart<<" rhs.astart: "<< rhs.astart<<"\n";
        return lhs.astart > rhs.astart;
        });



//for(int i =0 ;i < fData.uniBlocks.size();++i){
//std::cout<<fData.uniBlocks[i].achr<<"\n";
//}
    std::stable_sort(fData.uniBlocks.begin(), fData.uniBlocks.end(), [](const BLOCK& lhs, const BLOCK& rhs) -> bool{
  //  std::cout<<lhs.achr.compare(rhs.achr)<<"\n";
    if(lhs.achr.compare(rhs.achr) <= 0) return true;
//
//    std::cout<<"lhs.astart: " <<lhs.astart<<"  rhs.astart: " << rhs.astart<<"\n";
//    std::cout<<"lhs.achr: "<< lhs.achr<<"\n";
//    std::cout<<"rhs.achr: "<<rhs.achr<<"\n";

   // if(i ==0) return true;
//std::cout<<"done\n";
        return false;        });



    std::stable_sort(fData.dupBlocks.begin(), fData.dupBlocks.end(), [](const BLOCK &lhs, const BLOCK &rhs){
        return lhs.astart > rhs.astart ;
        });

    std::stable_sort(fData.dupBlocks.begin(), fData.dupBlocks.end(), [](const BLOCK &lhs, const BLOCK &rhs) -> bool{
        if(lhs.achr.compare(rhs.achr) <= 0) return true;
        return false;
        });


//    for(int i = 1000; i < 2000;++i){
//        std::cout<<i<<"  "<<fData.uniBlocks[i].astart<<"  "<<fData.uniBlocks[i].aend<<"  "<<fData.uniBlocks[i].achr<<"\n";
//    }

//     for(int i = 1000; i < 2000;++i){
//        std::cout<<i<<"  "<<fData.dupBlocks[i].astart<<"  "<<fData.dupBlocks[i].aend<<"  "<<fData.dupBlocks[i].achr<<"\n";
//    }

    std::cout<<"uSIZE: "<<fData.uniBlocks.size()<<"  mSize: "<<fData.dupBlocks.size()<<"\n";

	std::cout<<"FINISHED PARSE DUP"<<"\n";
	return(fData);
}


void filterBlocks(std::vector<BLOCK> &blocks, int threshold){
	// Identify duplicated regions in both genomes.
	// Ties are solved by comparing the chromosomes of the aligment.
	// If one block is CTX and other is not, then CTX block is
	// assigned as duplicated one. If both blocks are CTX or non-CTX
	// then length*identity score is used, with lower score block
	// being termed as duplicate. If the score is same, then the
	// later block is termed as duplicate.

	int size = blocks.size();

	// Duplicates in Genome A
	for(int i=0; i<size; i++){
		if (blocks[i].astate == DUP) continue;
		int check =0;
		for(int j=i+1; j< size and check==0; j++){
			if(blocks[i].achr.compare(blocks[j].achr)==0){
				if(blocks[i].astart > blocks[j].astart){
					printf("Error in input file. Positions in the first genome need to be sorted and larger \n"
							"i.astart: %d \t j.astart: %d \n", blocks[i].astart, blocks[j].astart);
					exit(1);
				}
				if(blocks[i].astart == blocks[j].astart){
					if(blocks[i].aend > blocks[j].aend){
						blocks[j].astate = DUP;
					}
					else{
						if(blocks[i].aend < blocks[j].aend){
							blocks[i].astate = DUP;
							check = 1;
						}
						else{
							if(blocks[i].iden > blocks[j].iden){
								blocks[j].astate = DUP;
							}
							else{
								blocks[i].astate = DUP;
								check = 1;
							}
						}
					}
				}
				else{
					if(blocks[i].aend > blocks[j].astart){
						if(blocks[i].aend >= blocks[j].aend){
							blocks[j].astate = DUP;
						}
						else{
							if(abs(blocks[i].astart - blocks[j].astart) < threshold and abs(blocks[i].aend - blocks[j].aend) < threshold){
								if(blocks[i].achr.compare(blocks[i].bchr)==0 and blocks[j].achr.compare(blocks[j].bchr)!=0){
									blocks[j].astate = DUP;
									continue;
								}
								if(blocks[i].achr.compare(blocks[i].bchr)!=0 and blocks[j].achr.compare(blocks[j].bchr)==0){
									blocks[i].astate = DUP;
									check = 1;
									continue;
								}
								float iScore = blocks[i].alen * blocks[i].iden;
								float jScore = blocks[j].alen * blocks[j].iden;
								if(iScore > jScore){
									blocks[j].astate = DUP;
								}
								else if(iScore < jScore){
									blocks[i].astate = DUP;
									check = 1;
								}
								else {
									blocks[j].astate = DUP;
									std::cerr << "Indistinguishable regions on Genome A:\n\t Chr: "<<blocks[i].achr<<"\t Block1 start: "<<blocks[i].astart<<
											"\t Block1 end: "<<blocks[i].aend<<"\t Block2 start: "<<blocks[j].astart<<"\t Block2 end: "<<blocks[j].aend<<"\n";
								}
							}else{
								if(abs(blocks[i].astart - blocks[j].astart) >threshold and abs(blocks[i].aend - blocks[j].aend) < threshold){
									blocks[j].astate = DUP;
								}
								else{
									if(abs(blocks[i].astart - blocks[j].astart) <threshold and abs(blocks[i].aend - blocks[j].aend) > threshold){
										blocks[i].astate = DUP;
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
	for(int i=0; i< size; i++){
		if(blocks[i].bstate == DUP) continue;
		int check = 0;
		for(int j = i+1; j < size and check == 0; j++){
			if(blocks[i].bchr.compare(blocks[j].bchr) == 0){
				if(blocks[i].bstart > blocks[j].bstart){
					if(blocks[i].bstart >= blocks[j].bend) continue;
					if (blocks[i].bend <= blocks[j].bend){
						blocks[i].bstate = DUP;
						check = 1;
					}
					else{
						if(abs(blocks[i].bstart - blocks[j].bstart) <threshold&& abs(blocks[i].bend - blocks[j].bend) < threshold){
							if(blocks[i].achr.compare(blocks[i].bchr)==0 and blocks[j].achr.compare(blocks[j].bchr)!=0){
								blocks[j].bstate = DUP;
								continue;
							}
							if(blocks[i].achr.compare(blocks[i].bchr)!=0 and blocks[j].achr.compare(blocks[j].bchr)==0){
								blocks[i].bstate = DUP;
								check = 1;
								continue;
							}
							float iScore = blocks[i].blen * blocks[i].iden;
							float jScore = blocks[j].blen * blocks[j].iden;

							if(iScore > jScore){
								blocks[j].bstate = DUP;
							}
							else if(iScore < jScore){
								blocks[i].bstate = DUP;
								check = 1;
							}
							else {
								blocks[j].bstate = DUP;
								std::cerr << "Indistinguishable regions on Genome B:\n\t Chr: "<<blocks[i].bchr<<"\t Block1 start: "<<blocks[i].bstart<<
										"\t Block1 end: "<<blocks[i].bend<<"\t Block2 start: "<<blocks[j].bstart<<"\t Block2 end: "<<blocks[j].bend<<"\n";
							}
						}
						else{
							if(abs(blocks[i].bstart - blocks[j].bstart) >threshold&& abs(blocks[i].bend - blocks[j].bend) < threshold){
								blocks[i].bstate = DUP;
								check = 1;
							}
							else{
								if(abs(blocks[i].bstart - blocks[j].bstart) <threshold&& abs(blocks[i].bend - blocks[j].bend) > threshold){
									blocks[j].bstate = DUP;
								}
							}
						}
					}
				}
				else{
					if(blocks[i].bstart == blocks[j].bstart){
						if(blocks[i].bend > blocks[j].bend){
							blocks[j].bstate = DUP;
						}
						else{
							if(blocks[i].bend < blocks[j].bend){
								blocks[i].bstate = DUP;
								check = 1;
							}
							else{
								if(blocks[i].iden > blocks[j].iden){
									blocks[j].bstate = DUP;
								}
								else{
									blocks[i].bstate = DUP;
									check = 1;
								}
							}
						}
					}
					else{
						if(blocks[i].bend > blocks[j].bstart){
							if(blocks[i].bend >= blocks[j].bend){
								blocks[j].bstate = DUP;
							}
							else{
								if(abs(blocks[i].bstart - blocks[j].bstart) <threshold&& abs(blocks[i].bend - blocks[j].bend) < threshold){
									if(blocks[i].achr.compare(blocks[i].bchr)==0 and blocks[j].achr.compare(blocks[j].bchr)!=0){
										blocks[j].bstate = DUP;
										continue;
									}
									if(blocks[i].achr.compare(blocks[i].bchr)!=0 and blocks[j].achr.compare(blocks[j].bchr)==0){
										blocks[i].bstate = DUP;
										check = 1;
										continue;
									}
									float iScore = blocks[i].blen * blocks[i].iden;
									float jScore = blocks[j].blen * blocks[j].iden;
									if(iScore > jScore){
										blocks[j].bstate = DUP;
									}
									else if(iScore < jScore){
										blocks[i].bstate = DUP;
										check = 1;
									}
									else {

										std::cout<<"iscore: "<<iScore<<"  jScore: "<<jScore<<"\n";

										blocks[j].bstate = DUP;
										std::cerr << "Indistinguishable regions on Genome B:\n\t Chr: "<<blocks[i].bchr<<"\t Block1 start: "<<blocks[i].bstart<<
												"\t Block1 end: "<<blocks[i].bend<<"\t Block2 start: "<<blocks[j].bstart<<"\t Block2 end: "<<blocks[j].bend<<"\n";
									}
								}
								else{
									if(abs(blocks[i].bstart - blocks[j].bstart) >threshold && abs(blocks[i].bend - blocks[j].bend) < threshold){
										blocks[j].bstate = DUP;
									}
									else{
										if(abs(blocks[i].bstart - blocks[j].bstart) <threshold && abs(blocks[i].bend - blocks[j].bend) > threshold){
											blocks[i].bstate = DUP;
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

	for (int i=0; i<size; i++){
		if(!(blocks[i].astate == DUP or blocks[i].bstate == DUP)){
			blocks[i].state = UNI;
			continue;
		}

		if((blocks[i].astate == DUP) xor (blocks[i].bstate == DUP)){
			blocks[i].state = DUP;
			continue;
		}

		if(blocks[i].astate == DUP and blocks[i].bstate == DUP){
			blocks[i].state = RED;
			continue;
		}


	}
	std::cout<<"FINISHED UNIQUE IDENTIFICATION \n";
}

void filterBlocks(std::vector<BLOCK> &mBlocks, std::vector<BLOCK> const &blocks, int threshold){

	int size = blocks.size();
	int msize = mBlocks.size();
int countA = 0;
int countB = 0;

	// Duplicates in Genome A
	for(int i=0; i<msize; ++i){
		for(int j=0; j< size; ++j){
			if(mBlocks[i].achr.compare(blocks[j].achr)==0 and blocks[j].state == UNI){
				if(!(((mBlocks[i].astart - blocks[j].astart) < -threshold ) or ((mBlocks[i].aend - blocks[j].aend) > threshold ))){
					mBlocks[i].astate = DUP;
					countA++;
					break;
				}
			}
		}
	}



	/*
				if(mBlocks[i].astart > blocks[j].astart){
					if(mBlocks[i].astart >= blocks[j].aend) continue;
					if (mBlocks[i].aend <= blocks[j].aend){
						mBlocks[i].astate = DUP;
						check = 1;
					}
					else{
						//if(!(abs(mBlocks[i].astart - blocks[j].astart) > 50 and abs(mBlocks[i].bend - blocks[j].bend) > 50)){
						if(!((mBlocks[i].bend - blocks[j].bend) > 50)){
							mBlocks[i].astate = DUP;
							check = 1;
						}
						else{
							if(abs(blocks[i].bstart - blocks[j].bstart) > 50 && abs(blocks[i].bend - blocks[j].bend) < 50){
								blocks[i].astate = DUP;
								check = 1;
							}
						}
					}
				}
				else{

				if(blocks[i].astart > blocks[j].astart){
					printf("Error in input file. Positions in the first genome need to be sorted and larger \n"
							"i.astart: %d \t j.astart: %d \n", blocks[i].astart, blocks[j].astart);
					exit(1);
				}
					if(mBlocks[i].astart == blocks[j].astart){
						if(!((mBlocks[i].aend - blocks[j].aend) > 50)){
				//			continue;
					//	}
				//		else{
							mBlocks[i].astate = DUP;
							if(mBlocks[i].aend < blocks[j].aend){
								mBlocks[i].astate = DUP;
								check = 1;
							}
							else{
								if(blocks[i].iden > blocks[j].iden){
									blocks[j].astate = DUP;
								}
								else{
									blocks[i].astate = DUP;
									check = 1;
								}
							}
						}
					}
					else{
						if(!(((mBlocks[i].astart - blocks[j].astart) < -50) or ((mBlocks[i].aend - blocks[j].aend) > 50))){
							if(mBlocks[i].aend >= blocks[j].aend){
								blocks[j].astate = DUP;
							}
							else{
								if(abs(blocks[i].astart - blocks[j].astart) < 50 && abs(blocks[i].aend - blocks[j].aend) < 50){
									if(blocks[i].state != DUP) blocks[i].astate = TEMP_DUP;
									if(blocks[j].state != DUP) blocks[j].astate = TEMP_DUP;
								}
								else{
									if(abs(blocks[i].astart - blocks[j].astart) > 50 && abs(blocks[i].aend - blocks[j].aend) < 50){
										blocks[j].astate = DUP;
									}
									else{
										if(abs(blocks[i].astart - blocks[j].astart) < 50 && abs(blocks[i].aend - blocks[j].aend) > 50){
											blocks[i].astate = DUP;
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
	 */
	//	Duplicates in Genome B

	for(int i=0; i< msize; ++i){
		for(int j = 0;j< size; ++j){
			if((mBlocks[i].bchr.compare(blocks[j].bchr) == 0) and blocks[j].state == UNI){
				if(!(((mBlocks[i].bstart - blocks[j].bstart) < -threshold ) or ((mBlocks[i].bend - blocks[j].bend) > threshold ))){
					mBlocks[i].bstate = DUP;
					countB++;
					break;
				}
			}
		}
	}
	for (int i=0; i<msize; i++){
		if(!(mBlocks[i].astate == DUP or mBlocks[i].bstate == DUP)){
			mBlocks[i].state = UNI;
			/*if(count>0){
				printf("i: %d   b.start: %d\n",i, mBlocks[i].bstart);

				count--;
			}
*/			continue;
		}

		if((mBlocks[i].astate == DUP) xor (mBlocks[i].bstate == DUP)){ // or (blocks[i].astate == NA and blocks[i].bstate == DUP)){
			mBlocks[i].state = DUP;
			continue;
		}

		if(mBlocks[i].astate == DUP and mBlocks[i].bstate == DUP){
			mBlocks[i].state = RED;
			continue;
		}
	}
	std::cout<<"FINISHED UNIQUE IDENTIFICATION FOR M_FILE \n";
}


/*

int compareABlocks(std::vector<BLOCK> &blocks, int i, int j){

		if(blocks[i].astate == DUP){
		compareABlocks(blocks, blocks[i].aDupOf, j);
	}
	if(blocks[j].astate == DUP){
		compareABlocks(blocks,i, blocks[j].bDupOf);
	}

	if(strcmp(blocks[i].achr, blocks[j].achr)==0) { // and blocks[j].state == UNI){
		if(blocks[i].astart > blocks[j].astart){
			if(blocks[i].astart >= blocks[j].aend) continue;
			if (blocks[i].aend <= blocks[j].aend){
				blocks[i].astate = DUP;
				check = 1;
			}
			else{
				if(abs(blocks[i].astart - blocks[j].astart) < 50 && abs(blocks[i].bend - blocks[j].bend) < 50){
					blocks[i].astate = TEMP_DUP;
				}
				else{
					if(abs(blocks[i].bstart - blocks[j].bstart) > 50 && abs(blocks[i].bend - blocks[j].bend) < 50){
						blocks[i].astate = DUP;
						check = 1;
					}
				}
			}
		}
		else{

					if(blocks[i].astart > blocks[j].astart){
						printf("Error in input file. Positions in the first genome need to be sorted and larger \n"
								"i.astart: %d \t j.astart: %d \n", blocks[i].astart, blocks[j].astart);
						exit(1);
					}
			if(blocks[i].astart == blocks[j].astart){
				if(blocks[i].aend > blocks[j].aend){
					blocks[j].astate = DUP;
				}
				else{
					if(blocks[i].aend < blocks[j].aend){
						blocks[i].astate = DUP;
						check = 1;
					}
					else{
						if(blocks[i].iden > blocks[j].iden){
							blocks[j].astate = DUP;
						}
						else{
							blocks[i].astate = DUP;
							check = 1;
						}
					}
				}
			}
			else{
				if(blocks[i].aend > blocks[j].astart){
					if(blocks[i].aend >= blocks[j].aend){
						blocks[j].astate = DUP;
					}
					else{
						if(abs(blocks[i].astart - blocks[j].astart) < 50 && abs(blocks[i].aend - blocks[j].aend) < 50){
							if(blocks[i].state != DUP) blocks[i].astate = TEMP_DUP;
							if(blocks[j].state != DUP) blocks[j].astate = TEMP_DUP;
						}
						else{
							if(abs(blocks[i].astart - blocks[j].astart) > 50 && abs(blocks[i].aend - blocks[j].aend) < 50){
								blocks[j].astate = DUP;
							}
							else{
								if(abs(blocks[i].astart - blocks[j].astart) < 50 && abs(blocks[i].aend - blocks[j].aend) > 50){
									blocks[i].astate = DUP;
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





for(int i = 0; i<BLOCK_NUM; i++){
	int check = 0;
	for(int j = i+1; j<BLOCK_NUM and check==0; j++){
		if(blocks[i].astart > blocks[j].astart and blocks[i].aend < blocks[j].aend and blocks[i].achr == blocks[j].achr){
			blocks[i].astate = DUP;
			check =1;
		}
	}

	check =0;
	for(int j = i+1; j <BLOCK_NUM and check==0; check++){
		if(blocks[i].bstart > blocks[j].bstate and blocks[i].bend < blocks[j].bend and blocks[i].bchr == blocks[j].achr){
			blocks[i].bstate = DUP;
			check =1;
		}
	}

	if(blocks[i].astate != DUP and blocks[i].bstate!=DUP){
		blocks[i].state = UNI;
	}
}

return(0);
}


 */


