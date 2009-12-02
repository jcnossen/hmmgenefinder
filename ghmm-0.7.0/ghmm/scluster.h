/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/scluster.h
*       Authors:  Bernhard Knab, Benjamin Georgi
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
*       This file is version $Revision: 1191 $ 
*                       from $Date: 2005-06-21 11:56:12 +0200 (Tue, 21 Jun 2005) $
*             last change by $Author: cic99 $.
*
*******************************************************************************/
#ifndef SCLUSTER_H
#define SCLUSTER_H
#include <ghmm/sequence.h>
#include <ghmm/smodel.h>
#include <ghmm/sreestimate.h>


#ifdef __cplusplus
extern "C" {
#endif

/**
   @name scluster
 */

/*@{
 */

#define CLASSIFY 0              /* Switch for Classificator: 0 == MD, 1 == MAW */

/**
   Cluster structure: All models and sequences. */
  struct scluster_t {
  /** 
  Vector of SHMMs pointers */
    smodel **smo;
  /** 
    Vector of sequence_t pointers; to store the sequences, that 
    belong to the models */
    sequence_d_t **smo_seq;
  /** 
    Number of models to read in */
    int smo_number;
  /** 
    Number of sequences for each model */
    long *seq_counter;
  /** 
    log(P) for the model */
    double *smo_Z_MD;
  /** a posteriori probability for the Modelle to calculate the objective
      fucntion in case of a MAW classificator. Is calculated using smap_bayes */
    double *smo_Z_MAW;
  };
/**
 */
  typedef struct scluster_t scluster_t;


/**
   Frees the memory allocated for a scluster_t struct.
   @return 0 for success; -1 for error
   @param scl pointer to scl struct
*/
  int scluster_t_free (scluster_t * scl);


/**
   Writes out the final models.
   @return 0 for success; -1 for error
   @param cl cluster of models to write
   @param sqd
   @param outfile output file
   @param out_filename name of the output file
 */
  int scluster_out (scluster_t * cl, sequence_d_t * sqd, FILE * outfile,
                    char *argv[]);

/**
   Avoids empty models going out as outputs by assigning them a random 
   sequence. This may lead to a produce of a new empty model - therefore
   change out sequences until a non empty model is found. (Quit after 100 
   iterations to avoid an infinite loop). 
   @return 0 for success; -1 for error
   @param sqd sequences to generate the model from
   @param cl cluster for the models
 */
  int scluster_avoid_empty_smodel (sequence_d_t * sqd, scluster_t * cl);

/**
   Makes a cluster and reestimates the HMMs.
   @return 0 for success; -1 for error
   @param argv vector of input files, one with sequences, one with models, 
   one for output and one with labels for the sequences - in this order.
 */
  int scluster_hmm (char *argv[]);

/**
   Updates the cluster with additional sequences.
   @return 0 for success; -1 for error
   @param cl cluster to update
   @param sqd sequences to update the cluster with
 */
  int scluster_update (scluster_t * cl, sequence_d_t * sqd);

/**
   Updates a label
   @return number of changes made
   @param oldlabel label to update
   @param up to date label for comparison
   @param seq_number number of sequences
   @param smo_changed tells, which labels have been changed
 */
  long scluster_update_label (long *oldlabel, long *seq_label,
                              long seq_number, long *smo_changed);

/**
   Prints out the likelihood values for the cluster.
   @param outfile output file
   @param cl cluster of models and sequences
 */
  void scluster_print_likelihood (FILE * outfile, scluster_t * cl);

/** 
    Determines form an already calculated probability matrix, which model 
    fits best to a certain sequence. 
    @return index of best model if success, otherwize -1
    @param cl cluster 
    @param seq_id ID of the sequence in question
    @param all_log_p matrix containing the probability of each sequence
    for each model
    @param log_p the probability of the sequence in question for the 
    best fitting model
*/
  int scluster_best_model (scluster_t * cl, long seq_id, double **all_log_p,
                           double *log_p);

/**
   Calculates the logarithmic probability of sequences for a model.
   @param cs sequences and model 
 */
  void scluster_prob (smosqd_t * cs);

/* int scluster_labels_from_kmeans(sequence_d_t *sqd, int smo_number); */

/**
   Creates random labels for a vector of sequences
   @return 0 for success; -1 for error
   @param sqd vector of sequences
   @param smo_number number of models (needed to determine the interval
   for the random numbers)
*/
  int scluster_random_labels (sequence_d_t * sqd, int smo_number);

/** Calculates the aposteriori prob. $\log(p(\lambda_best | O[seq\_id]))$, 
    where $\lambda_best$ is the model with highest apost. prob.
    @return 0 for success; -1 for error
    @param cl cluster
    @param sqd the sequence in question
    @param seq_id the ID of the sequence
    @param lob_apo the results
*/
  int scluster_log_aposteriori (scluster_t * cl, sequence_d_t * sqd,
                                int seq_id, double *log_apo);

/**
   Prints the input vector for scluster_hmm
 */
  void scluster_print_header (FILE * file, char *argv[]);

#ifdef __cplusplus
}
#endif
/*@} */
#endif                          /* SCLUSTER_H */
