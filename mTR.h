//
//  mTR.h
//  
//
//  Created by Shinichi Morishita on 2017/09/20.
//
//

// Key parameters

#define MAX_NUM_READS    100000   // The maximum number of all reads in a given fasta file. Each entry need 1KB.
#define MAX_INPUT_LENGTH 200000  // The maximum length of each read


//  Arrays for storing results
#define MIN_MATCH_RATIO 0.7      // The minimum threshold of match ratio
// between the estimated repeat unit and the repeat in a given raw read
#define BLK 1024                // Block size of input buffer.
#define MAX_PERIOD 500          // Maximum period length
#define MIN_REP_LEN_RATIO 0.6   // The threshold of the ratio of the actual repeat unit length to the estimated.
#define MIN_REP_LEN 20
#define MH_distance_threshold 0.2  // Two 2mer frequency distributions are identical if their Manhattan distance is less than or equal to this threshold.  A small threshold generates smaller groups of repeat units.


typedef struct {        // MAX_ID_LENGTH + MAX_EPRIOD + 28*4 = 612 bytes
    int     ID;  // 0,1,2,...
    char    readID[BLK];
    int     inputLen;
    int     rep_start;
    int     rep_end;
    int     repeat_len;
    int     rep_period;
    int     Num_freq_unit;
    int     Num_matches;
    int     Num_mismatches;
    int     Num_insertions;
    int     Num_deletions;
    int     Kmer;
    int     ConsensusMethod;     // 0 = progressive multiple alignment, 1 = De Bruijn graph search
    char    string[MAX_PERIOD];
    int     predicted_rep_period;
    int     freq_2mer[16];
} repeat_in_read;

repeat_in_read *repeats_in_all_reads;

// For clustering repeats

#define MIN_NUM_repTR 1

typedef struct TR1{    // Sorted according to actual_rep_period, freq_2mer, and Num_freq_unit
    int     ID;     //  20 B
    int     rep_period;
    int     freq_2mer[16];
    int     Num_freq_unit;
    int     repID;
    struct TR1 *repTR;      // This allows us to refer the struct itself.
    int     freq;
} TR;

TR      *TR_list;
TR   *repTR_list;

// External functions
#define MAX(a, b) ((a) > (b) ? (a) : (b))

int feed_rr_into_repeats_in_all_reads(char *inputFile);
void print_one_repeat_in_read(repeat_in_read rr);
void k_means_clustering(int read_cnt);
