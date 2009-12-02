/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/cluster.h
*       Authors:  Bernd Wichern
*
*       Copyright (C) 1998-2004 Alexander Schliep
*       Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
*       Copyright (C) 2002-2004 Max-Planck-Institut fuer Molekulare Genetik,
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

#ifndef CLUSTER_H
#define CLUSTER_H
#include <ghmm/sequence.h>
#include <ghmm/model.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 @name cluster
 */

/*@{ */

/**
 */
  struct cluster_t {
  /** 
    Vector of HMM-Model pointers */
    model **mo;
  /** 
    Vector of sequence_t pointers: for saving the sequence data
    that belongs to the models
  */
    sequence_t **mo_seq;
  /** 
    Number of models read in */
    int mo_number;
  };
/**
*/
  typedef struct cluster_t cluster_t;

/**
   Writess out the final models.
   @return 0 for success; -1 for error
   @param cl cluster of models to write
   @param sq 
   @param outfile output file
   @param out_filename name of the output file
 */
  int cluster_out (cluster_t * cl, sequence_t * sq, FILE * outfile,
                   char *out_filename);

/**
   Prevents that empty models are sent out (no associated seqences) by 
   associating a random sequence. Since it's possible to produce an empty model
   in this way, the sequences are shifted until a nonempty model is produced. (This 
   could be a never-ending process and therefore it's only done 100 times).
   @return 0 for success; -1 for error
   @param seq_label vector of labels for the sequences
   @param seq_number number of sequences
   @param mo_number number of models
 */
  int cluster_avoid_empty_model (long *seq_label, long seq_number,
                                 int mo_number);

/**
   Makes a cluster and reestimates the HMMs.
   @return 0 for success; -1 for error
   @param seq_file file of sequences
   @param mo_file file of initial models
   @param out_file output file
 */
  int cluster_hmm (char *seq_file, char *mo_file, char *out_file);

/**
   Updates the cluster with additional sequences.
   @return 0 for success; -1 for error
   @param cl cluster to update
   @param sq sequences to update the cluster with
 */
  int cluster_update (cluster_t * cl, sequence_t * sq);

/**
   Updates a label
   @return number of changes made
   @param oldlabel label to update
   @param up to date label for comparison
   @param seq_number number of sequences
 */
  long cluster_update_label (long *oldlabel, long *seq_label,
                             long seq_number);

/**
   Prints out the likelihood values for the cluster.
   @param outfile output file
   @param cl cluster of models and sequences
 */
  void cluster_print_likelihood (FILE * outfile, cluster_t * cl);

/*@} cluster documentation */

#ifdef __cplusplus
}
#endif
#endif                          /* CLUSTER_H */
