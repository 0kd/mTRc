//
//  k_means_clustering.c
//  
//
//  Created by Shinichi Morishita on 2017/10/06.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mTR.h"

int cmp_TR(TR aTR, TR bTR, int mode){
    // mode = 0: Ordered by <rep_period, freq_2mer> of aTR and bTR
    // mode = 1: Ordered by <rep_period, freq_2mer, Num_freq_unit> of aTR and bTR
    // mode = 2: Ordered by <frequency, identifier> of the representative TRs of aTR and bTR
    
    int diff;
    if(mode == 0 || mode == 1){
        // Compute the differenece between the unit lengths (rep_period)
        diff = aTR.rep_period - bTR.rep_period;
        // If the lengths of the two units are equal, compute the difference between the 2mer distribution (freq_2mer)
        if(diff != 0){
            return(diff);   // negative: <, positive means: >
        }else{ // diff == 0
            for(int i=0; i<16; i++){
                diff = aTR.freq_2mer[i] - bTR.freq_2mer[i];
                if(diff != 0){
                    return(diff);   // negative: <, positive means: >
                }
            }
            // If aTR.freq_2mer == bTR.freq_2mer; namely, aTR.freq_2mer[i] == bTR.freq_2mer[i] for all i.
            if( mode == 1){
                diff = aTR.Num_freq_unit - bTR.Num_freq_unit;
                return(diff);
            }else{
                return(0);
            }
        }
    }else if(mode == 2){
        // Ordering according to the frequency and identifier of the representative TR
        diff = -(aTR.repTR->freq - bTR.repTR->freq);
        if(diff == 0){
            diff = aTR.repTR->ID - bTR.repTR->ID;
            return(diff);
        }else{
            return(diff);
        }
    }else{
        return(0);
    }
}

int partition_TR_list(int left, int right, int mode) {
    // Select a number between left and right at random.
    int random = rand() % (right-left+1) + left;
    // Exchange target[right] and target[random].
    TR tmp = TR_list[right]; TR_list[right] = TR_list[random]; TR_list[random] = tmp;
    
    TR pivot = TR_list[right];
    int i = left-1; // i scans the array from the left.
    int j = right;  // j scans the array from the right.
    for (;;) {
        // Move from the left until hitting a value no less than the pivot.
        for(i++; cmp_TR(TR_list[i], pivot, mode) < 0; i++){}
        // Move from the right until hitting a value no greater than the pivot.
        for(j--; cmp_TR(pivot, TR_list[j], mode) < 0 && i < j; j--){}
        if (i >= j)  break;
        // Exchange target[i] and target[j].
        tmp = TR_list[i];  TR_list[i] = TR_list[j];  TR_list[j] = tmp;
    }
    // Exchange target[i] and target[right].
    tmp = TR_list[i];  TR_list[i] = TR_list[right];  TR_list[right] = tmp;
    return i;
}

// We set "mode" to 1 in the first call of this function.
// Afterwards, we set "mode" to 2 in the second call.
void sort_TR_List(int left, int right, int mode){
    if (left < right) {
        int i = partition_TR_list(left, right, mode); // i: Position of the pivot.
        sort_TR_List(left, i - 1, mode);
        sort_TR_List(i + 1, right, mode);
    }
}

int select_repTR_list(int Num_qualified_TRs){
    int freq = 1;
    int j = 0;
    int i;
    for(i=0; i<(Num_qualified_TRs - 1); i++){
        if(cmp_TR(TR_list[i], TR_list[i+1], 0) == 0){
            // mode = 0: Ordered by <rep_period, freq_2mer>
            // If two consecutive TRs are identical in terms of <rep_period, freq_2mer>,
            // put them into an identical group, and increment its size, denoted by freq.
            freq++;
        }else if(MIN_NUM_repTR <= freq){    // End of a block of identical TRs
            repTR_list[j] = TR_list[i];    // Set the last TR in one group to the reresentative
            repTR_list[j].freq = freq;
            for(int k=i; i-freq+1<=k; k--){   // Update representative IDs
                TR_list[k].repID =  repTR_list[j].ID;
                TR_list[k].repTR = &repTR_list[j];
            }
            j++;
            freq = 1;
        }
    }
    if(i == Num_qualified_TRs - 1 && MIN_NUM_repTR <= freq){
        repTR_list[j] = TR_list[i];
        repTR_list[j].freq = freq;
        for(int k=i; i-freq+1<=k; k--){   // Update representative IDs
            TR_list[k].repID =  repTR_list[j].ID;
            TR_list[k].repTR = &repTR_list[j];
        }
        j++;
    }
    return(j); // Number of representative TRs
}

int TRs_in_neighborhood(TR aTR, TR repTR){
    // return +1 if they are in neighborhood, and -1 otherwise
    int diff = 0;
    for(int i=0; i<16; i++){
        diff += abs(aTR.freq_2mer[i] - repTR.freq_2mer[i]);
    }
    // We here set the threshold to a temporary value, but it should be revised.
    if( MH_distance_threshold * repTR.rep_period < diff)
        return(-1);
    else
        return(1);
}

void revise_repTR_list(int Num_repTRs){
    // Revise the clustering to merge groups whose 2mer distributions are similar according to the measure defined in TRs_in_neighborhood. Select the largest group as the representative.
    for(int i=0; i<Num_repTRs; i++){
        TR aTR = repTR_list[i];
        // Search the list of TRs for those similar to the focal TR.
        // Candidate TRs are of length within 10% of the focal TR length.
        int lb = aTR.rep_period - (int)(aTR.rep_period * 0.1);
        int ub = aTR.rep_period + (int)(aTR.rep_period * 0.1);
        int max_freq = aTR.freq;
        int max_i = i;
        
        // Search backward
        for(int j=i-1; 0<=j; j--){  // Search positions
            TR repTR = repTR_list[j];
            if(lb <= repTR.rep_period){
                if(TRs_in_neighborhood(aTR, repTR) == 1 && max_freq < repTR.freq){
                    max_freq = repTR.freq;
                    max_i = j;
                }
            }else{
                break;
            }
        }
        // Search forward
        for(int j=i+1; j<Num_repTRs; j++){
            TR repTR = repTR_list[j];
            if(repTR.rep_period <= ub){
                if(TRs_in_neighborhood(aTR, repTR) == 1 && max_freq < repTR.freq){
                    max_freq = repTR.freq;
                    max_i = j;
                }
            }else{
                break;
            }
        }
        repTR_list[i].repID =  repTR_list[max_i].ID;
        repTR_list[i].repTR = &repTR_list[max_i];
    }
    
    for(int i=0; i<Num_repTRs; i++){
        TR *aTR = &repTR_list[i];
        if(aTR->ID != aTR->repID){
            // Search the root of representative groups recursively.
            while(aTR->ID != aTR->repID){
                aTR = aTR->repTR;
            }
            repTR_list[i].repID = aTR->ID;
            repTR_list[i].repTR = aTR;
            aTR->freq += repTR_list[i].freq;
        }
    }
}

void update_repTRs_in_TR_list(int Num_TRs){
    for(int i=0; i<Num_TRs; i++){
        TR *aTR = TR_list[i].repTR;
        if(aTR->ID != aTR->repID){
            // Search the root of representative groups recursively.
            while(aTR->ID != aTR->repID){
                aTR = aTR->repTR;
            }
            TR_list[i].repID = aTR->ID;
            TR_list[i].freq  = aTR->freq;
            TR_list[i].repTR = aTR;
        }
    }
}


void print_one_TR(TR aTR){
    
    char message[2000] = "";
    
    if(aTR.rep_period > 0){
        if(aTR.repTR == NULL){
            sprintf(message, "%i,\tUnit len = %i,\t#Units = %i,\trepID = %i,\tfreq = %i,\t2mers = ",
                    aTR.ID,
                    aTR.rep_period,
                    aTR.Num_freq_unit,
                    aTR.repID,
                    aTR.freq
                    );
        }else{
            sprintf(message, "%i,\tUnit len = %i,\t#Units = %i,\trepID = %i,\tfreq = %i,\t2mers = ",
                    aTR.ID,
                    aTR.rep_period,
                    aTR.Num_freq_unit,
                    aTR.repTR->ID,
                    aTR.repTR->freq
                    );
        }
        char tmp_msg[6]="";
        for(int i=0; i<16; i++){
            sprintf(tmp_msg, "%i,", aTR.freq_2mer[i]);
            strcat(message, tmp_msg);
        }
        strcat(message, "\tstring = ");
        strcat(message, repeats_in_all_reads[aTR.ID].string);
    }
    
    printf("%s\n", message);
    
}

void k_means_clustering(int read_cnt){
    
    TR_list      = (TR*) malloc(read_cnt*sizeof(TR));
    if( TR_list == NULL ){
        fprintf(stderr, "cannot allocate space for one of global variables in the heap.\n");
        exit(EXIT_FAILURE);
    }    
    repTR_list   = (TR*) malloc(read_cnt*sizeof(TR));
    if( repTR_list == NULL ){
        fprintf(stderr, "cannot allocate space for one of global variables in the heap.\n");
        exit(EXIT_FAILURE);
    }

    
    repeat_in_read aRR;
    
    // Select qualified TRs
    int j = 0;
    for(int i = 0; i < read_cnt; i++){
        aRR = repeats_in_all_reads[i];
        
        float match_ratio = (float)aRR.Num_matches / aRR.repeat_len;
        
        if(MIN_REP_LEN < (aRR.rep_period * aRR.Num_freq_unit)  &&
           MIN_MATCH_RATIO < match_ratio &&
           1 < aRR.Num_freq_unit )
        { // qualifed
            TR_list[j].ID = aRR.ID;
            TR_list[j].repID = aRR.ID;
            TR_list[j].repTR = NULL;
            TR_list[j].freq = 1;
            
            TR_list[j].rep_period = aRR.rep_period;
            for(int k =0; k<16; k++){
                TR_list[j].freq_2mer[k] = aRR.freq_2mer[k];
            }
            TR_list[j].Num_freq_unit = aRR.Num_freq_unit;
            // initialization

            j++;
        }
    }
    int Num_qualified_TRs = j;
    
    // Sort according to rep_period, freq_2mer, and Num_freq_unit by setting mode to 1
    sort_TR_List(0, Num_qualified_TRs - 1, 1);
    
    // Cluster identical TRs in terms of rep_period and freq_2mer
    // Determine the size of each group
    int Num_repTRs = select_repTR_list(Num_qualified_TRs);
    revise_repTR_list(Num_repTRs);

    update_repTRs_in_TR_list(Num_qualified_TRs);
    
    // Sort according to frequency and repID by setting mode to 2
    sort_TR_List(0, Num_qualified_TRs - 1, 2);
    
    int aRepID = 0;
    
    for(int j = 0; j < Num_qualified_TRs; j++){
        repeat_in_read tmp_rr = repeats_in_all_reads[ TR_list[j].ID ];
     
        if(aRepID != TR_list[j].repID){
            aRepID = TR_list[j].repID;
            printf("#members = %d\tunit-length = %d",
                   TR_list[j].repTR->freq,
                   TR_list[j].repTR->rep_period);
         
            //printf("\t2mer-frequency-vector = ");
            //for(int i=0; i<16; i++){ printf("%d ", TR_list[j].freq_2mer[i]);}
         
            printf("\n");
        }
        
        print_one_repeat_in_read( tmp_rr );

        if(j % BLK == 0){
            fflush( stdout );
        }
    }

    free(TR_list);
    free(repTR_list);
    
}

