# mTRc
clustering tandem repeats

Usage:  mTRc <list of tandem repeats>

<list of tandem repeats>
    char    *readID;
    int     inputLen;         The length of an input read 
    int     rep_start;        The start position of a predicted tandem repeat (TR) in the read
    int     rep_end;          The end position of a predicted TR in the read
    int     repeat_len;       The length of a predicted TR
    int     rep_period;       The length of the unit of a predicted TR
    int     Num_freq_unit;    The frequency of the unit of a predicted TR
    int     Num_matches;      The number of matches between the predicted TR and its corresponding region in the read 
    int     Num_mismatches;   The number of mismatches
    int     Num_insertions;   The number of insertions
    int     Num_deletions;    The number of deletions
    int     Kmer;             The length of k-mer used for predicting the TR
    int     ConsensusMethod;  
    char    *string;          The string of the repeat unit
  
 Output groups of tandem repeats
 Each group starts with a line of the form:
 #members = <number of TRs in the group>  unit-length = <length of the unit of the representative TR in the group>
 A list of tandem repeats in the group follows.
  
