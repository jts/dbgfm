#include <stdio.h>
#include "fm_index.h"

int main()
{
    printf("starting main\n");
    std::string test_bwt = "/home/jsimpson/simpsonlab/data/references/human.chr20.bwt";
    FMIndex index(test_bwt);
}
