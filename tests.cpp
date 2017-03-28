#include "tests.h"

void estimateThreshold(std::vector<BLOCK> &blocks){

    std::ofstream fout;
    fout.open("uniqueCount.tsv");

   // std::vector<int> uniCount;
    int size = blocks.size();

    for(int i = 0; i<1000;){

    int count =0;

    std::vector<BLOCK> temp_blocks = blocks;


    filterBlocks(temp_blocks, i);

    for(int j = 0; j < size; ++j){
    if(temp_blocks[j].state == UNI) count++;
    }

   // uniCount.push_back(count);
   fout<<i<<"\t"<<count<<"\n";
   i+=100;
    }




fout.close();
}
