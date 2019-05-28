# mTRc

clustering tandem repeats

Usage:  mTRc [list of tandem repeats in input reads]

One tandem repeat (TR) in a read has the form:

    char    *readID;
    int     inputLen;         The length of an input read 
    int     rep_start;        The start position of a predicted TR in the read
    int     rep_end;          The end position of a predicted TR in the read
    int     repeat_len;       The length of a predicted TR
    int     rep_period;       The length of the unit of a predicted TR
    int     Num_freq_unit;    The frequency of the unit of a predicted TR
    int     Num_matches;      The number of matches between the predicted TR and its corresponding region in the read 
    int     Num_mismatches;   The number of mismatches
    int     Num_insertions;   The number of insertions
    int     Num_deletions;    The number of deletions
    int     Kmer;             The length of k-mer used for predicting the TR
    int     ConsensusMethod;  1 = De Bruijn graph search, 0 = progressive multiple alignment
    char    *string;          The string of the repeat unit

Example:

    read1   5084    81  1559    1479    26  57  1239    0.837728    195     45      64      5       1   TGACTCTGGCCGTTCACCAAATTTAG   
  
The program outputs groups of highly similar tandem repeats. 

    Each group starts with a line of the form:
        #members = <number of TRs in the group>  unit-length = <length of the unit of the representative TR in the group>
    A list of tandem repeats in the group follows.
    Example:
    #members = 406  unit-length = 26
    read1   5084    81  1559    1479    26  57  1239    0.837728    195     45      64      5       1   TGACTCTGGCCGTTCACCAAATTTAG      
    read2   9191    81  1524    1444    26  58  1184    0.819945    231     29      116     6       1   TGACTCTGGCCGTTCACCAAATTTAG      
    read3   2007    91  2004    1914    26  75  1731    0.904389    159     24      69      5       1   TGACTCTGGCCGTTCACCAAATTTAG  
