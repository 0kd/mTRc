//---------------------------------------------------------------------
// Finding tandem repeats in long noizy reads
// Initial codes are developed by Shinichi Morishita and ....
//---------------------------------------------------------------------

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mTR.h"

void print_error_message(){
    fprintf(stderr, "mTRc <reads with repeats>: Cluster reads with repeats. Feed <reads with repeats> and print <clustering results>. \n");
}

int main(int argc, char *argv[])
{
    // mTRc <reads with repeats>
    //                              Feed <reads with repeats> and print <clustering results>.
    
    char *inputFile;
    int pretty_print = 0;
    
    if(argc == 2){  // Valid arguments
    }else{
        print_error_message();
        exit(EXIT_FAILURE);
    }

    //  Allocate space for global variables in the heap
    repeats_in_all_reads = (repeat_in_read*) malloc(MAX_NUM_READS*sizeof(repeat_in_read));
    if(repeats_in_all_reads == NULL){
        fprintf(stderr, "Fatal error: cannot allocate a space for repeats_in_all_reads\n");
        exit(EXIT_FAILURE);
    }
    
    // Cluster reads with repeats into similar groups
    int read_cnt = feed_rr_into_repeats_in_all_reads(argv[1]);
    //for(int i=0; i<read_cnt; i++){print_one_repeat_in_read( repeats_in_all_reads[i] );}
        
    k_means_clustering(read_cnt);
    
    
    //  Free space for global variables in the heap
    free(repeats_in_all_reads);

    return EXIT_SUCCESS;
}
