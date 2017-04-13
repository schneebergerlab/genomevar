/*
 * getDuplicates.cpp
 *
 *  Created on: Mar 16, 2017
 *      Author: goel
 */
#include "getDuplicates.h"

void parseDUP(std::vector<BLOCK> &blocks,std::vector<BLOCK> &mBlocks, int threshold, FILTEREDDATA &fData){
	std::cout<< "FINDING DUPLICATES\n";

	int usize = blocks.size();
    int msize = mBlocks.size();

	filterBlocks(blocks, threshold);

	for(int i = 0; i < usize; ++i){
		if(blocks[i].state == UNI){
			fData.uniBlocks.push_back(blocks[i]);
        }
        else if(blocks[i].state == DUP){
            fData.dupBlocks.push_back(blocks[i]);
        }
    }

	filterBlocks(mBlocks, fData.uniBlocks, threshold);

	std::vector<BLOCK> mUni;

	for(int i =0; i < msize;++i){
        if(mBlocks[i].state == UNI){
			mUni.push_back(mBlocks[i]);
			continue;
		}
        if(mBlocks[i].state == DUP) fData.dupBlocks.push_back(mBlocks[i]);
	}

	int mUniSize = mUni.size();

    filtermUni(mUni, threshold);
    //filterBlocks(mUni, threshold);


    int countUni = 0;
    int countDup = 0;
    int countRed = 0;

    for(int i = 0; i < mUniSize; ++i){
        if(mUni[i].state == UNI) ++countUni;
        if(mUni[i].state == DUP) ++countDup;
        if(mUni[i].state == RED) ++countRed;


    }

    std::cout<<"UNIQUE MUNI : "<<countUni<<"  "<<countDup<<"  "<<countRed<<"\n";


  //  annotateDup(fData.dupBlocks, fData.uniBlocks, threshold);
    for(int i = 0; i < mUniSize; ++i){

        if(mUni[i].state == UNI){
     //       std::cout<<mUni[i].astart<<" "<<mUni[i].aend<<" "<<mUni[i].achr <<" "<<mUni[i].bstart<<" "<<mUni[i].bend<<" "<<mUni[i].bchr<<"\n";
            fData.uniBlocks.push_back(mUni[i]);
        }
        else if (mUni[i].state == DUP){
            fData.dupBlocks.push_back(mUni[i]);
        }
    }

    sortBlocks(fData.uniBlocks);
    sortBlocks(fData.dupBlocks);

//    std::stable_sort(fData.uniBlocks.begin(), fData.uniBlocks.end(), []( const BLOCK &lhs, const BLOCK &rhs){
//        return lhs.astart > rhs.astart;
//        });
//
//    std::stable_sort(fData.uniBlocks.begin(), fData.uniBlocks.end(), [](const BLOCK& lhs, const BLOCK& rhs) -> bool{
//    if(lhs.achr.compare(rhs.achr) <= 0) return true;
//        return false;        });
//
//    std::stable_sort(fData.dupBlocks.begin(), fData.dupBlocks.end(), [](const BLOCK &lhs, const BLOCK &rhs){
//        return lhs.astart > rhs.astart ;
//        });
//    std::stable_sort(fData.dupBlocks.begin(), fData.dupBlocks.end(), [](const BLOCK &lhs, const BLOCK &rhs) -> bool{
//        if(lhs.achr.compare(rhs.achr) <= 0) return true;
//        return false;
//        });


    std::cout<<"uSIZE: "<<fData.uniBlocks.size()<<"  mSize: "<<fData.dupBlocks.size()<<"\n";

//	            std::cout<<"Size: "<<fData.uniBlocks.size()<<"\n";
//
//	            	filterBlocks(fData.dupBlocks, fData.dupBlocks, threshold);
//
//    int a =0;
//    int b =0;
//    int c=0;
//for(int i = 0; i< fData.dupBlocks.size();++i){
//    if(fData.dupBlocks[i].state == UNI){
//        a++;
//    }else if(fData.dupBlocks[i].state == DUP){
//        b++;
//    }else c++;
//}
//
//std::cout<<a<<"  "<<b<<"  "<<c<<"\n";

	annotateDup(fData.dupBlocks, fData.uniBlocks, threshold);


    std::ofstream mout;
    mout.open("uniBlocks.txt");
    for(int i = 0; i < fData.uniBlocks.size(); ++i){
        mout<<fData.uniBlocks[i].astart<<"\t"<<fData.uniBlocks[i].aend<<"\t"<<fData.uniBlocks[i].bstart<<"\t"<<
            fData.uniBlocks[i].bend<<"\t"<<fData.uniBlocks[i].dir<<"\t"<<fData.uniBlocks[i].iden<<"\t"<<
            fData.uniBlocks[i].alen<<"\t"<<fData.uniBlocks[i].blen<<"\t"<<fData.uniBlocks[i].achr<<"\t"<<
            fData.uniBlocks[i].bchr<<"\t"<<fData.uniBlocks[i].astate<<"\t"<<fData.uniBlocks[i].bstate<<"\t"<<
            fData.uniBlocks[i].state<<"\n";
    }
    mout.close();

	std::vector<BLOCK> filterTest = mblocks;


	for(int i =0; i < filterTest.size(); ++i){
	    filterTest[i].state = NA;
	    filterTest[i].state = NA;
	    filterTest[i].state = NA;
	}
	int countAlarm = 0;
	filtermUni(filterTest, threshold);

    std::vector<int> countData = getStateCount(filterTest);

    std::cout<<"STATE COUNT :  "<<countData[0] <<"  "<<countData[1]<<"  "<<countData[2]<<"\n";




  //  std::ofstream mout;
    mout.open("mBlocks.txt");
    for(int i = 0; i < filterTest.size(); ++i){
        mout<<filterTest[i].astart<<"\t"<<filterTest[i].aend<<"\t"<<filterTest[i].bstart<<"\t"<<
            filterTest[i].bend<<"\t"<<filterTest[i].dir<<"\t"<<filterTest[i].iden<<"\t"<<
            filterTest[i].alen<<"\t"<<filterTest[i].blen<<"\t"<<filterTest[i].achr<<"\t"<<
            filterTest[i].bchr<<"\t"<<filterTest[i].astate<<"\t"<<filterTest[i].bstate<<"\t"<<
            filterTest[i].state<<"\n";
    }
    mout.close();

	std::vector<BLOCK> testUni, testDup;

	for(int i = 0; i < filterTest.size(); ++i){
        if(filterTest[i].state == UNI){
            testUni.push_back(filterTest[i]);
            continue;

        }
        if(filterTest[i].state == DUP){
            testDup.push_back(filterTest[i]);
            continue;
        }
    }

    filterBlocks(testUni, threshold);
    countData = getStateCount(testUni);
    std::cout<<"STATE COUNT :  "<<countData[0] <<"  "<<countData[1]<<"  "<<countData[2]<<"\n";
  //  annotateDup(fData.dupBlocks, fData.uniBlocks, threshold);

    sortBlocks(testUni);
    sortBlocks(testDup);

    annotateDup(testDup, testUni, threshold);
    fData.uniBlocks.clear();
    fData.uniBlocks.insert(fData.uniBlocks.begin(), testUni.begin(), testUni.end());
    fData.dupBlocks.clear();
    fData.dupBlocks.insert(fData.dupBlocks.begin(), testDup.begin(), testDup.end());

    printDup(testDup, testUni);

	std::cout<<"FINISHED PARSE DUP"<<"\n";
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
	int countA = 0, countB = 0, countC = 0;

	// Duplicates in Genome A
	for(int i=0; i<size; i++){
		//if (blocks[i].astate == DUP) continue;
		if (blocks[i].state == UNI) continue;
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
//		if(blocks[i].bstate == DUP) continue;
		if(blocks[i].astart == 9438785 and blocks[i].bstart == 10641459 and blocks[i].achr == "Chr4"){
            std::cout<<"CHECKING\n";
        }
		int check = 0;
		for(int j = i+1; j < size and check == 0; j++){
			if(blocks[i].bchr.compare(blocks[j].bchr) == 0){
				if(blocks[i].bstart > blocks[j].bstart){
					if(blocks[i].bstart >= blocks[j].bend) continue;
					if (blocks[i].bend <= blocks[j].bend){
						blocks[i].bstate = DUP;
						if(blocks[i].astart == 9438785 and blocks[i].bstart == 10641459 and blocks[i].achr == "Chr4"){
							std::cout<<blocks[j].astart<<" "<<blocks[j].aend<<" "<<blocks[j].achr <<" "<<blocks[j].bstart<<" "<<blocks[j].bend<<" "<<blocks[j].bchr<<"\n";
                            std::cout<<"FOUND 1 \n";
							        }

						check = 1;
					}
					else{
						if(abs(blocks[i].bstart - blocks[j].bstart) < threshold and abs(blocks[i].bend - blocks[j].bend) < threshold){
							if(blocks[i].achr.compare(blocks[i].bchr)==0 and blocks[j].achr.compare(blocks[j].bchr)!=0){
								blocks[j].bstate = DUP;
								continue;
							}
							if(blocks[i].achr.compare(blocks[i].bchr)!=0 and blocks[j].achr.compare(blocks[j].bchr)==0){
								blocks[i].bstate = DUP;
								if(blocks[i].astart == 9438785 and blocks[i].bstart == 10641459 and blocks[i].achr == "Chr4"){
						//	std::cout<<blocks[j].astart<<" "<<blocks[j].aend<<" "<<blocks[j].achr <<" "<<blocks[j].bstart<<" "<<blocks[j].bend<<" "<<blocks[j].bchr<<"\n";
                            std::cout<<"FOUND 2 \n";
							        }
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
								if(blocks[i].astart == 9438785 and blocks[i].bstart == 10641459 and blocks[i].achr == "Chr4"){
						//	std::cout<<blocks[j].astart<<" "<<blocks[j].aend<<" "<<blocks[j].achr <<" "<<blocks[j].bstart<<" "<<blocks[j].bend<<" "<<blocks[j].bchr<<"\n";
                            std::cout<<"FOUND 3 \n";
							        }
								check = 1;
							}
							else {
								blocks[j].bstate = DUP;
								std::cerr << "Indistinguishable regions on Genome B:\n\t Chr: "<<blocks[i].bchr<<"\t Block1 start: "<<blocks[i].bstart<<
										"\t Block1 end: "<<blocks[i].bend<<"\t Block2 start: "<<blocks[j].bstart<<"\t Block2 end: "<<blocks[j].bend<<"\n";
							}
						}
						else{
							if(abs(blocks[i].bstart - blocks[j].bstart) > threshold and abs(blocks[i].bend - blocks[j].bend) < threshold){
								blocks[i].bstate = DUP;
								if(blocks[i].astart == 9438785 and blocks[i].bstart == 10641459 and blocks[i].achr == "Chr4"){
						//	std::cout<<blocks[j].astart<<" "<<blocks[j].aend<<" "<<blocks[j].achr <<" "<<blocks[j].bstart<<" "<<blocks[j].bend<<" "<<blocks[j].bchr<<"\n";
                            std::cout<<"FOUND 4 \n";
							        }
								check = 1;
							}
							else{
								if(abs(blocks[i].bstart - blocks[j].bstart) < threshold and abs(blocks[i].bend - blocks[j].bend) > threshold){
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
								if(blocks[i].astart == 9438785 and blocks[i].bstart == 10641459 and blocks[i].achr == "Chr4"){
						//	std::cout<<blocks[j].astart<<" "<<blocks[j].aend<<" "<<blocks[j].achr <<" "<<blocks[j].bstart<<" "<<blocks[j].bend<<" "<<blocks[j].bchr<<"\n";
                            std::cout<<"FOUND 5 \n";
							        }

								check = 1;
							}
							else{
								if(blocks[i].iden > blocks[j].iden){
									blocks[j].bstate = DUP;
								}
								else{
									blocks[i].bstate = DUP;
									if(blocks[i].astart == 9438785 and blocks[i].bstart == 10641459 and blocks[i].achr == "Chr4"){
						//	std::cout<<blocks[j].astart<<" "<<blocks[j].aend<<" "<<blocks[j].achr <<" "<<blocks[j].bstart<<" "<<blocks[j].bend<<" "<<blocks[j].bchr<<"\n";
                            std::cout<<"FOUND 6 \n";
							        }

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
								if(abs(blocks[i].bstart - blocks[j].bstart) <threshold and abs(blocks[i].bend - blocks[j].bend) < threshold){
									if(blocks[i].achr.compare(blocks[i].bchr)==0 and blocks[j].achr.compare(blocks[j].bchr)!=0){
										blocks[j].bstate = DUP;
										continue;
									}
									if(blocks[i].achr.compare(blocks[i].bchr)!=0 and blocks[j].achr.compare(blocks[j].bchr)==0){
										blocks[i].bstate = DUP;
										if(blocks[i].astart == 9438785 and blocks[i].bstart == 10641459 and blocks[i].achr == "Chr4"){
						//	std::cout<<blocks[j].astart<<" "<<blocks[j].aend<<" "<<blocks[j].achr <<" "<<blocks[j].bstart<<" "<<blocks[j].bend<<" "<<blocks[j].bchr<<"\n";
                            std::cout<<"FOUND 7 \n";
							        }

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
										if(blocks[i].astart == 9438785 and blocks[i].bstart == 10641459 and blocks[i].achr == "Chr4"){
						//	std::cout<<blocks[j].astart<<" "<<blocks[j].aend<<" "<<blocks[j].achr <<" "<<blocks[j].bstart<<" "<<blocks[j].bend<<" "<<blocks[j].bchr<<"\n";
                            std::cout<<"FOUND 8 \n";
							        }

										check = 1;
									}
									else {

									//	std::cout<<"iscore: "<<iScore<<"  jScore: "<<jScore<<"\n";

										blocks[j].bstate = DUP;
										std::cerr << "Indistinguishable regions on Genome B:\n\t Chr: "<<blocks[i].bchr<<"\t Block1 start: "<<blocks[i].bstart<<
												"\t Block1 end: "<<blocks[i].bend<<"\t Block2 start: "<<blocks[j].bstart<<"\t Block2 end: "<<blocks[j].bend<<"\n";
									}
								}
								else{
									if(abs(blocks[i].bstart - blocks[j].bstart) >threshold and abs(blocks[i].bend - blocks[j].bend) < threshold){
										blocks[j].bstate = DUP;
									}
									else{
										if(abs(blocks[i].bstart - blocks[j].bstart) <threshold and abs(blocks[i].bend - blocks[j].bend) > threshold){
											blocks[i].bstate = DUP;
											if(blocks[i].astart == 9438785 and blocks[i].bstart == 10641459 and blocks[i].achr == "Chr4"){
						//	std::cout<<blocks[j].astart<<" "<<blocks[j].aend<<" "<<blocks[j].achr <<" "<<blocks[j].bstart<<" "<<blocks[j].bend<<" "<<blocks[j].bchr<<"\n";
                            std::cout<<"FOUND 9 \n";
							        }

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
			++countA;
			continue;
		}

		if((blocks[i].astate == DUP) xor (blocks[i].bstate == DUP)){
			blocks[i].state = DUP;
			++countB;
			continue;
		}

		if(blocks[i].astate == DUP and blocks[i].bstate == DUP){
			blocks[i].state = RED;
			++countC;
			continue;
		}
	}
	std::cout<<size<<"  "<<countA<< "  "<<countB<<"  "<<countC<<"   FINISHED UNIQUE IDENTIFICATION \n";
}

void filterBlocks(std::vector<BLOCK> &mBlocks, std::vector<BLOCK> const &blocks, int threshold){
	int usize = blocks.size();
	int msize = mBlocks.size();
    int countA = 0;
    int countB = 0;
    int countC = 0;

	// Duplicates in Genome A
	for(int i=0; i<msize; ++i){
		for(int j=0; j< usize; ++j){
			if(mBlocks[i].achr.compare(blocks[j].achr)==0 and blocks[j].state == UNI){
				if(!(((mBlocks[i].astart - blocks[j].astart) < -threshold ) or ((mBlocks[i].aend - blocks[j].aend) > threshold ))){
					mBlocks[i].astate = DUP;
					//countA++;
					break;
				}
			}
		}
	}

//	Duplicates in Genome B

	for(int i=0; i< msize; ++i){
		for(int j = 0;j< usize; ++j){
			if((mBlocks[i].bchr.compare(blocks[j].bchr) == 0) and blocks[j].state == UNI){
				if(!(((mBlocks[i].bstart - blocks[j].bstart) < -threshold ) or ((mBlocks[i].bend - blocks[j].bend) > threshold ))){
					mBlocks[i].bstate = DUP;
					//countB++;
					break;
				}
			}
		}
	}
	for (int i=0; i<msize; ++i){
//	    if(mBlocks[i].astart == 26590117 and mBlocks[i].aend == 26592092 and mBlocks[i].achr == "Chr5" and mBlocks[i].bstart == 25985707){
//	        std::cout<<"state  "<<mBlocks[i].astate<<"  "<<mBlocks[i].bstate<<"\n";
//	    }
		if(!(mBlocks[i].astate == DUP or mBlocks[i].bstate == DUP)){
			mBlocks[i].state = UNI;
			++countA;
			continue;
		}

		if((mBlocks[i].astate == DUP) xor (mBlocks[i].bstate == DUP)){
			mBlocks[i].state = DUP;
			++countB;
			continue;
		}

		if(mBlocks[i].astate == DUP and mBlocks[i].bstate == DUP){
			mBlocks[i].state = RED;
			++countC;
			continue;
		}
	}
	std::cout<<msize<<"  "<<countA<< "  "<<countB<<"  "<<countC<<"   FINISHED UNIQUE IDENTIFICATION FOR M_FILE \n";
}

void filtermUni(std::vector<BLOCK> &blocks, int threshold){

    int msize = blocks.size();
    int naCount = 0;
    int uniCount = 0;

    std::vector<BLOCK> uniBlocks, dupBlocks, redBlocks;
    std::vector<BLOCK> temp;

    for(int i =0; i < msize ;++i){
        blocks[i].astate = NA;
        blocks[i].bstate = NA;
        blocks[i].state = NA;
    }

    do{
        naCount = 0;
        uniCount = 0;

        filterBlocks(blocks, threshold);

        for(int j = 0; j < blocks.size(); ++j){
            if(blocks[j].state == UNI){
                uniBlocks.push_back(blocks[j]);
          //      blocks.erase(blocks.begin() + j);
        //        --j;
                ++uniCount;
                continue;
            }
            blocks[j].astate = NA;
            blocks[j].bstate = NA;
            blocks[j].state = NA;
            temp.push_back(blocks[j]);
        }

        blocks.clear();
        blocks.insert(blocks.end(), temp.begin(), temp.end());
        temp.clear();




        if(uniCount == 0){
            std::cout<<"Unresolvable alignments found. Saved in unresolved.txt\n";
            std::ofstream unr;
            unr.open("unresolved.txt");
            for(int j = 0; j < blocks.size(); ++j){
                unr<<blocks[j].astart<<"\t"<<blocks[j].aend<<"\t"<<blocks[j].bstart<<"\t"<<
                    blocks[j].bend<<"\t"<<blocks[j].dir<<"\t"<<blocks[j].iden<<"\t"<<
                    blocks[j].alen<<"\t"<<blocks[j].blen<<"\t"<<blocks[j].achr<<"\t"<<
                    blocks[j].bchr<<"\n";
            }

            blocks.clear();
            unr.close();
            break;
        }

        sortBlocks(uniBlocks);
//        std::stable_sort(uniBlocks.begin(), uniBlocks.end(), []( const BLOCK &lhs, const BLOCK &rhs){
//            return lhs.astart > rhs.astart;
//        });
//
//        std::stable_sort(uniBlocks.begin(), uniBlocks.end(), [](const BLOCK& lhs, const BLOCK& rhs) -> bool{
//            if(lhs.achr.compare(rhs.achr) <= 0) return true;
//            return false;
//        });

        filterBlocks(blocks, uniBlocks, threshold);
        for(int j = 0; j < blocks.size(); ++j){
            if(blocks[j].state == DUP){
                dupBlocks.push_back(blocks[j]);
                //blocks.erase(blocks.begin() + j);
                //--j;
                continue;
            }
            else if(blocks[j].state == RED){
                redBlocks.push_back(blocks[j]);
                //blocks.erase(blocks.begin() + j);
                //--j;
                continue;
            }else if(blocks[j].state == UNI){
                blocks[j].state = NA;
                temp.push_back(blocks[j]);
                ++naCount;
            }
        }

        blocks.clear();
        blocks.insert(blocks.end(), temp.begin(), temp.end());
        temp.clear();
//        for(int i = 0; i < msize; ++i){
//            if(blocks[i].state == NA) ++naCount;
//        }
    }while(naCount != 0);


/*
    for(int i = 0; i < msize; ++i){
        for(int j=i+1; j < msize; j++){
            if(abs(blocks[i].astart - blocks[j].astart) + abs(blocks[i].aend - blocks[j].aend) < threshold){
                blocks[i].astate = DUP;
                blocks[j].astate = DUP;
                continue;
            }
            if(blocks[i].astart <= blocks[j].astart  and  blocks[i].aend >= blocks[j].aend){
                    blocks[j].astate = DUP;
            }
            else if(blocks[i].astart >= blocks[j].astart  and  blocks[i].aend <= blocks[j].aend){
                blocks[i].astate = DUP;
            }
        }
    }


    for(int i = 0; i < msize; ++i){
        for(int j=i+1; j < msize; j++){
            if(blocks[i].bchr.compare(blocks[j].bchr)==0){
                if(abs(blocks[i].bstart - blocks[j].bstart) + abs(blocks[i].bend - blocks[j].bend) < threshold){
                    blocks[i].bstate = DUP;
                    blocks[j].bstate = DUP;
                    continue;
                }
                if(blocks[i].bstart <= blocks[j].bstart  and  blocks[i].bend >= blocks[j].bend){
                        blocks[j].bstate = DUP;
                }
                else if(blocks[i].bstart >= blocks[j].bstart  and  blocks[i].bend <= blocks[j].bend){
                    blocks[i].bstate = DUP;
                }
            }
        }
    }

    std::vector<int> uniIndices;

    for (int i=0; i<msize; ++i){
		if(!(blocks[i].astate == DUP or blocks[i].bstate == DUP)){
			blocks[i].state = UNI;
						uniIndices.push_back(i);

			continue;
		}
		else{
            blocks[i].astate = 0;
            blocks[i].bstate = 0;
            continue;
		}
    }

    int change = 0;
    do{
            std::cout<<"**START***";
        change = 0;
        for(int i = 0; i < blocks.size(); ++i){

            if(blocks[i].state == NA){
                for(int j = 0; j < blocks.size(); ++j){
                    if(blocks[j].state == UNI){
                        if(blocks[i].achr.compare(blocks[j].achr)==0){
                            if(!(((blocks[i].astart - blocks[j].astart) < -threshold ) or ((blocks[i].aend - blocks[j].aend) > threshold ))){
                                blocks[i].astate = DUP;
                            }
                        }
                        if((blocks[i].bchr.compare(blocks[j].bchr) == 0)){
                            if(!(((blocks[i].bstart - blocks[j].bstart) < -threshold ) or ((blocks[i].bend - blocks[j].bend) > threshold ))){
                                blocks[i].bstate = DUP;
                            }
                        }
                    }
                    else if(blocks[j].state == DUP){
                        if(abs(blocks[i].astart - blocks[j].astart) + abs(blocks[i].aend - blocks[j].aend) < threshold){
                            blocks[i].astate = DUP;
                        }else if(blocks[i].astart >= blocks[j].astart and blocks[i].aend <= blocks[j].aend){
                            blocks[i].astate = DUP;
                        }

                        if(abs(blocks[i].bstart - blocks[j].bstart) + abs(blocks[i].bend - blocks[j].bend) < threshold){
                            blocks[i].bstate = DUP;
                        }else if(blocks[i].bstart >= blocks[j].bstart and blocks[i].bend <= blocks[j].bend){
                            blocks[i].bstate = DUP;
                        }
                    }
                }
            }
        }

        for (int i=0; i<blocks.size(); ++i){
            if(blocks[i].state == NA){
                if((blocks[i].astate == DUP) xor (blocks[i].bstate == DUP)){
                    blocks[i].state = DUP;
                    ++change;
                    continue;
                }

                if(blocks[i].astate == DUP and blocks[i].bstate == DUP){
                    blocks.erase(blocks.begin()+i);
                    --i;
                }
            }
        }

        std::cout<<"CHANGE: "<<change<<"\n";
    }while(change!=0);
*/
    std::cout<<"mUniSize: "<<blocks.size()<<"\n";


blocks.insert( blocks.end(), uniBlocks.begin(), uniBlocks.end() );
blocks.insert( blocks.end(), dupBlocks.begin(), dupBlocks.end() );
blocks.insert( blocks.end(), redBlocks.begin(), redBlocks.end() );
}


void annotateDup(std::vector<BLOCK> &dupBlocks, std::vector<BLOCK> &uniBlocks, int threshold){
std::cout<<"ANNOTATING\n";
    int msize = dupBlocks.size();
    int usize = uniBlocks.size();
    int cnt = 0;
    int inI = 0;
    int check;
    int missed = 0;


    for(int i = 0; i < msize; ++i){
		inI++;
		check  = 0;
        for(int j = 0; j < usize; ++j){

        	if(dupBlocks[i].achr.compare(uniBlocks[j].achr) == 0){

				if(!(((dupBlocks[i].astart - uniBlocks[j].astart) < -threshold ) or ((dupBlocks[i].aend - uniBlocks[j].aend) > threshold ))){
					dupBlocks[i].dupOf = j;
					uniBlocks[j].duplicates.push_back(i);
					cnt++;
					check = 1;
					break;
				}
			}
			if(dupBlocks[i].bchr.compare(uniBlocks[j].bchr) == 0){
				if(!(((dupBlocks[i].bstart - uniBlocks[j].bstart) < -threshold ) or ((dupBlocks[i].bend - uniBlocks[j].bend) > threshold ))){
					dupBlocks[i].dupOf = j;
					uniBlocks[j].duplicates.push_back(i);
					cnt++;
					check = 1;
					break;
				}
			}
        }

        if(check == 0){
        missed++;
    //    std::cout<<dupBlocks[i].astart<<" "<<dupBlocks[i].bstart<<" "<<dupBlocks[i].achr<<"\n";
        //std::cout<<i<<"\n";
        }
    }
    std::cout<<inI<<"  "<<cnt<<"  "<<missed<<"  ANNOTATING FINISHED\n";
}

void printDup(std::vector<BLOCK> const &dupBlocks, std::vector<BLOCK> const &uniBlocks){
    int usize = uniBlocks.size();
    int msize = dupBlocks.size();

    for(int i = 0; i < usize; ++i){
        if(uniBlocks[i].duplicates.size() != 0){
            dupOutFile<<"#DUP\t"<<uniBlocks[i].astart<<"\t"<<uniBlocks[i].aend<<"\t"<<uniBlocks[i].achr<<
                " - "<<uniBlocks[i].bstart<<"\t"<<uniBlocks[i].bend<<"\t"<<uniBlocks[i].bchr<<"\n";
            for(int j = 0; j < uniBlocks[i].duplicates.size(); ++j){
                dupOutFile<< dupBlocks[uniBlocks[i].duplicates[j]].astart<<"\t"<<dupBlocks[uniBlocks[i].duplicates[j]].aend<<
                    "\t"<<dupBlocks[uniBlocks[i].duplicates[j]].bstart<<"\t"<<dupBlocks[uniBlocks[i].duplicates[j]].bend<<
                    "\t"<<dupBlocks[uniBlocks[i].duplicates[j]].dir<<"\t"<<dupBlocks[uniBlocks[i].duplicates[j]].achr<<
                    "\t"<<dupBlocks[uniBlocks[i].duplicates[j]].bchr<<"\n";
            }
        }
    }
}


std::vector<int> getStateCount(std::vector<BLOCK> const &blocks){
    int bsize = blocks.size();
    int countUNI = 0;
    int countDUP = 0;
    int countRED = 0;

    for( int i = 0; i < bsize; ++i){
        if(blocks[i].state == UNI) ++countUNI;
        if(blocks[i].state == DUP) ++countDUP;
        if(blocks[i].state == RED) ++countRED;
    }
    std::vector<int> v = {countUNI, countDUP, countRED};
    return v;

}

void sortBlocks(std::vector<BLOCK>& blocks){
    std::stable_sort(blocks.begin(), blocks.end(), []( const BLOCK &lhs, const BLOCK &rhs){
        return lhs.astart > rhs.astart;
    });

    std::stable_sort(blocks.begin(), blocks.end(), [](const BLOCK& lhs, const BLOCK& rhs) -> bool{
        if(lhs.achr.compare(rhs.achr) <= 0) return true;
        return false;
    });
}
