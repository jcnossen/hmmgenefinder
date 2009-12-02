/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/sequence.h
*       Authors:  Bernd Wichern, Benjamin Georgi
*
*       Copyright (C) 1998-2004 Alexander Schliep 
*       Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
*	Copyright (C) 2002-2004 Max-Planck-Institut fuer Molekulare Genetik, 
*                               Berlin
*                                   
*       Contact: schliep@ghmm.org             
*
*       This library is free software; you can redistribute it and/or
*       modify it under the terms of the GNU Library General Public
*       License as published by the Free Software Foundation; either
*       version 2 of the License, or (at your option) any later version.
*
*       This library is distributed in the hope that it will be useful,
*       but WITHOUT ANY WARRANTY; without even the implied warranty of
*       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*       Library General Public License for more details.
*
*       You should have received a copy of the GNU Library General Public
*       License along with this library; if not, write to the Free
*       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*
*
*       This file is version $Revision: 1253 $ 
*                       from $Date: 2005-08-10 11:39:24 +0200 (Wed, 10 Aug 2005) $
*             last change by $Author: cic99 $.
*
*******************************************************************************/
#ifndef SEQUENCE_H
#define SEQUENCE_H

#ifdef __cplusplus
extern "C" {
#endif

/**@name sequences  (double and int) */
/*@{ (Doc++-Group: sequence) */
/** @name struct sequence_t
    Sequence structure for integer sequences. 
    Contains an array of sequences and corresponding
    data like sequence label, sequence weight, etc. Sequences may have different
    length.    
 */

  struct sequence_t {
  /** sequence array. sequence[i] [j] = j-th symbol of i-th seq.
   */
    int **seq;

  /** matrix of state ids, can be used to save the viterbi path during sequence generation.
   ATTENTION: is NOT allocated by sequence_calloc  */
    int **states;

  /** array of sequence length */
    int *seq_len;
  /**  array of sequence labels */
    long *seq_label;
  /**  array of sequence IDs*/
    double *seq_id;
  /** positiv! sequence weights.  default is 1 = no weight */
    double *seq_w;
  /** total number of sequences */
    long seq_number;
  /** sum of sequence weights */
    double total_w;

    /* matrix of state labels corresponding to seq */
    int **state_labels;
    /* number of labels for each sequence */
    int *state_labels_len;
  };
  typedef struct sequence_t sequence_t;

/** @name struct sequence_d_t
    Sequence structure for double sequences. 
    Contains an array of sequences and corresponding
    data like sequnce label, sequence weight, etc. Sequences may have different
    length.    
 */
  struct sequence_d_t {
  /** sequence array. sequence[i][j] = j-th symbol of i-th seq. */
    double **seq;
  /** array of sequence length */
    int *seq_len;
  /**  array of sequence labels */
    long *seq_label;
  /**  array of sequence IDs*/
    double *seq_id;
  /** positive! sequence weights.  default is 1 = no weight */
    double *seq_w;
  /** total number of sequences */
    long seq_number;
  /** sum of sequence weights */
    double total_w;
  };
  typedef struct sequence_d_t sequence_d_t;


#ifdef __cplusplus
}
#endif
/* don't include model.h at the beginning of this file. struct sequence_t has
   to be known in model.h */
#include <ghmm/model.h>
#include <ghmm/smodel.h>
#ifdef __cplusplus
extern "C" {
#endif

/**
   Reads one or several arrays of integer sequences. 
   Calls sequence\_read\_alloc, where reading
   and memory allocation is done. 
   @return pointer to sequence array
   @param filename    input filename
   @param seq\_arrays number of sequence arrays read
*/
  sequence_t **sequence_read (const char *filename, int *seq_arrays);

/**
   Reading of one integer sequence field. Memory alloc here.
   @param s scanner
   @return array of sequences
*/
  sequence_t *sequence_read_alloc (scanner_t * s);

/**
   Reads one or several arrays of double sequences. 
   Calls sequence\_read\_alloc, where reading
   and memory allocation is done. 
   @return pointer to sequence array
   @param filename    input filename
*/
  sequence_d_t **sequence_d_read (const char *filename, int *sqd_number);

/**
   Reading of one double sequence field. Memory alloc here.
   @param s scanner
   @return array of sequences
*/
  sequence_d_t *sequence_d_read_alloc (scanner_t * s);

/** Truncate double sequences in a given sequence array. 
    Useful for Testing;
   @return truncated sqd_field; 
   @param sqd\_in sequence arrays for truncation
   @param sqd\_arrays number of sequence arrays
   @param  trunc\_ratio 0 means  no truncation, 1 max. truncation
   @param seed rng seed
*/

  sequence_d_t **sequence_d_truncate (sequence_d_t ** sqd_in, int sqd_arrays,
                                      double trunc_ratio, int seed);



/**
  Extract a single sequence from a larger sequence_t into a new struct.
  
  @return sequence_t struct containing a single sequence
  @param sq   source sequence_t
  @param index   index of sequence to extract
*/
sequence_t *sequence_get_singlesequence(sequence_t *sq, int index);

/**
  Extract a single sequence_d from a larger sequence_d_t into a new struct.
  
  @return sequence_d_t struct containing a single sequence
  @param sq   source sequence_d_t
  @param index   index of sequence to extract
*/
sequence_d_t *sequence_d_get_singlesequence(sequence_d_t *sq, int index);

/*XXX TEST: frees everything but the seq field */

/**
  Free a sequence_t struct which holds as sequence a reference to a sequence in a different
  sequence_t. The function deallocates everything but the reference.
*/
int sequence_subseq_free (sequence_t ** sq);

/**
  Free a sequence_d_t struct which holds as sequence a reference to a sequence in a different
  sequence_d_t. The function deallocates everything but the reference.
*/
int sequence_d_subseq_free (sequence_d_t ** sqd);






/** Generates all possible integer sequence of lenght n from an alphabet with
    M letters. Use lexicographical ordering. Memory allocation here.
    @param n      length of sequences
    @param M     size of alphabet
    @return array of generated integer sequences
*/
  sequence_t *sequence_lexWords (int n, int M);

/**
   Determine best model for a given integer sequence. 
   Choose from the set of models the 
   one with the highest likelihood for the given sequence.
   @param mo            array of models
   @param model\_number  number of models
   @param sequence      sequence
   @param seq\_len      length of sequence
   @param log\_p         log likelihood of the sequence given the best model
   @return index of best\_model (between 0 and model\_number - 1)
*/
  int sequence_best_model (model ** mo, int model_number, int *sequence,
                           int seq_len, double *log_p);

/**
   Make sure that the sequences only contain allowed symbols. 
   (between 0 and max\_symbol - 1)
   @param sq          sequences
   @param max_symb    number of different symbols
   @return            -1 for error, 0 for no errors
*/
  int sequence_check (sequence_t * sq, int max_symb);

/**
  copy one integer sequence. Memory for target has to be allocated outside.
  @param target  target sequence
  @param source source sequence
  @param len     length of source sequence
  */
  void sequence_copy (int *target, int *source, int len);

/**
  copy one double sequence. Memory for target has to be allocated outside.
  @param target  target sequence
  @param source source sequence
  @param len     length of source sequence
  */
  void sequence_d_copy (double *target, double *source, int len);

/**
  Adds all integer sequences, sequence lengths etc 
  from source to target. Memory allocation is done here.
  @param target target sequence structure
  @param source  source sequence structure
  @return -1 for error, 0 for success
  */
  int sequence_add (sequence_t * target, sequence_t * source);


/**
  Adds all double sequences, sequence lengths etc 
  from source to target. Memory allocation is done here.
  @param target target sequence structure
  @param source  source sequence structure
  @return -1 for error, 0 for success
  */
  int sequence_d_add (sequence_d_t * target, sequence_d_t * source);

/**
  Prints one array of integer sequences in a file.
  @param file       output file
  @param sequence    array of sequences
  */
  void sequence_print (FILE * file, sequence_t * sequence);

/**
  Prints one array of integer sequences in a xml file
  @param file       output file
  @param sequence   array of sequences
  */
  void sequence_print_xml (FILE * file, sequence_t * sequence);

/**
   Prints one array of integer sequences in Mathematica format.
   (List of lists)
   @param file       output file
   @param sq    array of sequences
   @param name arbitrary sequence name for usage in Mathematica.
 */
  void sequence_mathematica_print (FILE * file, sequence_t * sq, char *name);

/**
  Prints one array of double sequences in a file.
  @param file       output file
  @param sqd    array of sequences
  @param discrete   switch: 0 means double output for symbols,  
     1 means truncate symbols to integer
  */
  void sequence_d_print (FILE * file, sequence_d_t * sqd, int discrete);

/**
   Prints one array of double sequences in Mathematica format.
   (List of lists)
   @param file       output file
   @param sqd    array of sequences
   @param name arbitrary sequence name for usage in Mathematica.
 */
  void sequence_d_mathematica_print (FILE * file, sequence_d_t * sqd,
                                     char *name);

/** Output of double sequences suitable for gnuplot. One symbol per line,
    sequences seperated by double newline.
    @param file output file
    @param sqd array of double sequences
*/
  void sequence_d_gnu_print (FILE * file, sequence_d_t * sqd);

/**
   Cleans integer sequence pointers in sequence struct. sets 
   seq\_number to zero.
   Differs from sequence\_free since memory is not freed here. 
   @param sq sequence structure
  */
  void sequence_clean (sequence_t * sq);

/**
   Cleans double sequence pointers in sequence struct. sets 
   seq\_number to zero.
   Differs from sequence\_free since memory is not freed here. 
   @param sqd sequence structure
  */
  void sequence_d_clean (sequence_d_t * sqd);

/**
  Frees all memory in a given array of integer sequences.
  @param sq sequence  structure
  @return 0 for succes, -1 for error
  */
  int sequence_free (sequence_t ** sq);

/**
  Frees all memory in a given array of double sequences.
  @param sq sequence  structure
  @return 0 for succes, -1 for error
  */
  int sequence_d_free (sequence_d_t ** sq);

/**
   Return biggest symbol in an interger sequence.
   @param sq sequence structure
   @return max value
 */
  int sequence_max_symbol (sequence_t * sq);

/**
   Memory allocation for an integer sequence struct. Allocates arrays of lenght
   seq\_number. NO allocation for the actual sequence, since its length is 
   unknown.
   @param seq\_number:  number of sequences
   @return:     pointer of sequence struct
*/
  sequence_t *sequence_calloc (long seq_number);

/**
   Memory allocation for a double  sequence struct. Allocates arrays of lenght
   seq\_number. NO allocation for the actual sequence, since its length is 
   unknown.
   @param seq\_number:  number of sequences
   @return:     pointer of sequence struct
*/
  sequence_d_t *sequence_d_calloc (long seq_number);

/**
   Copies array of integer sequences to double sequences.
   @return       double sequence struct (target)
   @param sq    integer sequence struct (source)
   */
  sequence_d_t *sequence_d_create_from_sq (const sequence_t * sq);

/**
   Copies array of double sequences into an array of integer
   sequences. Truncates positions after decimal point.
   @return       integer sequence struct (target)
   @param sq    double sequence struct (source)
   */
  sequence_t *sequence_create_from_sqd (const sequence_d_t * sqd);

/** 
    Determines max sequence length in a given int sequence struct.
    @author Peter Pipenbacher
    @param sqd sequence struct
    @return max sequence length
 */
  int sequence_max_len (const sequence_t * sqd);

/** 
    Determines max sequence length in a given double sequence struct.
    @param sqd sequence struct
    @return max sequence length
 */
  int sequence_d_max_len (const sequence_d_t * sqd);

/**
  Calculates a mean sequence of a given array of double sequences.
  Missing values of shorter sequences a assumed to be zero.
  @param sqd sequence struct
  @return pointer of sequence struct containing the mean sequence
  */
  sequence_d_t *sequence_d_mean (const sequence_d_t * sqd);

/**
   Calculates the scatter matrix of an array of double sequences. 
   Missing parts of short sequences are NOT taken into account.
   @return        scatter matrix
   @param sqd     sequence struct
   @param sqd     (calculated) dimension of scatter matrix
  */
  double **sequence_d_scatter_matrix (const sequence_d_t * sqd, int *dim);

/**
   Calculates transition class for a given double sequence
   at a specified position. Very application specific!!! Currently 
   implemented only dummy function: allways returns 0 which
   means no usage of multiple transition classes.
   @param O double sequence
   @param index position for class calculation
   @param osum sum of symbols upto index
   @return currently always 0
 */
  int sequence_d_class (const double *O, int index, double *osum);
  /*int sequence_d_class(const sequence_d_t *sqd, const int seq_number, int index, double *osum, );*/

/** Divides randomly a given array of double sequences into two sets. 
    Useful if a training and test set is needed. Memory allocation is done 
    here.
    @param sqd input sequence array
    @param sqd\_train training sequences
    @param sqd\_test test sequences
    @param train\_ratio ratio of number of train vs number of test sequences
    @return 0 for success, -1 for error
*/
  int sequence_d_partition (sequence_d_t * sqd, sequence_d_t * sqd_train,
                            sequence_d_t * sqd_test, double train_ratio);


/** 
    Copies all entries from one sequence in a source array to a target array.
    No memory allocation here.
    @param target double sequence target
    @param source double sequence source
    @param t_num position in target array
    @param s_num position in source array
*/
  void sequence_d_copy_all (sequence_d_t * target, long t_num,
                            sequence_d_t * source, long s_num);

/** Log-Likelihood function in a mixture model:
    (mathe mode?)
    $\sum_k w^k \log( \sum_c (\alpha_c p(O^k | \lambda_c)))$
    @param smo pointer to array of smodels
    @param smo\_number number of models
    @param sqd sequence struct
    @param like log likelihood
*/
  int sequence_d_mix_like (smodel ** smo, int smo_number, sequence_d_t * sqd,
                           double *like);

#ifdef __cplusplus
}
#endif
#endif
/*@} (Doc++-Group: sequence) */
