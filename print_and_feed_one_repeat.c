//
//  print_one_repeat.c
//  
//
//  Created by Shinichi Morishita on 2017/10/06.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mTR.h"

int char2int(char c){
    int val;
    switch(c){
        case 'A': val = 0; return(val);
        case 'C': val = 1; return(val);
        case 'G': val = 2; return(val);
        case 'T': val = 3; return(val);
        default: fprintf(stderr, "fatal input char %c", c); exit(EXIT_FAILURE);
    }
}

void freq_2mer_array(char* st, int len, int *freq_2mer){
    for(int i=0; i<16; i++){
        freq_2mer[i] = 0;
    }
    int index;
    for(int i=1; i<len; i++){
        index = char2int(st[i-1]) * 4 + char2int(st[i]);
        freq_2mer[index]++;
    }
    // wrap around and concatenate the last and first characters
    index = char2int(st[len-1]) * 4 + char2int(st[0]);
    freq_2mer[index]++;
}

int feed_rr_into_repeats_in_all_reads(char *inputFile){
    FILE *fp = fopen(inputFile, "r");
    if(fp == NULL){
        fprintf(stderr, "fatal error: cannot open %s\n", inputFile);
        exit(EXIT_FAILURE);
    }
    //fprintf(stderr, "Input = %s\n", inputFile);
    
    char tmp_s[BLK];
    int i = 0;
    repeat_in_read tmp;
    
    while( fscanf(fp, "%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%s\t%i\t%i\t%i\t%i\t%i\t%s\n",
             tmp.readID,
             &tmp.inputLen,
             &tmp.rep_start,
             &tmp.rep_end,
             &tmp.repeat_len,
             &tmp.rep_period,
             &tmp.Num_freq_unit,
             &tmp.Num_matches,
             tmp_s,
             &tmp.Num_mismatches,
             &tmp.Num_insertions,
             &tmp.Num_deletions,
             &tmp.Kmer,
             &tmp.ConsensusMethod,
             tmp.string) != EOF )
    {
        freq_2mer_array(tmp.string, tmp.rep_period, tmp.freq_2mer);
        repeats_in_all_reads[i] = tmp;
        repeats_in_all_reads[i].ID = i;
        i++;
    }
    
    //fprintf(stderr, "Number of input reads = %i\n", i);
    fclose(fp);
    return(i);
}

void print_one_repeat_in_read(repeat_in_read rr){
    char message[2000] = "";
    sprintf(message, "%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%f\t%i\t%i\t%i\t%i\t%i\t%s\t",
            rr.readID,
            rr.inputLen,
            rr.rep_start,
            rr.rep_end,
            rr.repeat_len,
            rr.rep_period,
            rr.Num_freq_unit,
            rr.Num_matches,
            (double)rr.Num_matches/rr.repeat_len,
            rr.Num_mismatches,
            rr.Num_insertions,
            rr.Num_deletions,
            rr.Kmer,
            rr.ConsensusMethod,
            rr.string
            );
    /*
    char tmp_msg[6]="";
    for(int i=0; i<16; i++){
        sprintf(tmp_msg, "%d ", rr.freq_2mer[i]);
        strcat(message, tmp_msg);
    }
     */
    printf("%s\n", message);
}


