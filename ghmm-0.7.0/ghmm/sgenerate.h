/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/sgenerate.h
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
#ifndef SGENERATE_H
#define SGENERATE_H


#include <ghmm/smodel.h>
#include <ghmm/sequence.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
   @name generation and extension of sequences from shmm
*/

/*@{
 */

/**
 Help "modus" for sgenerate_extensions.
*/

  typedef enum {
    viterbi_viterbi,
    viterbi_all,
    all_viterbi,
    all_all
  } sgeneration_mode_t;

/**
   Makes part sequences longer given a model. There are some different possibilities
   to do this (only Viterbi or take all paths into account, combined with 
   sequence begin and end).

   Which method to use:
   0 = viterbi\_viterbi, 
   1 = viterbi\_all, 
   2 = all\_viterbi, 
   3 = all\_all
   (zunaechst nur all\_all moeglich)
   @return pointer to a vector of all sequences (given initial sequence and generated 
   end sequence)
   @param smo:         given Model
   @param sqd_short:   vector of initial sequences
   @param seed:        initial value for random value generator (int)
   @param global_len:  wanted length of sequences (=0: automatically over final states)
   @param mode:        which method to use for the generator
 */
  sequence_d_t *sgenerate_extensions (smodel * smo, sequence_d_t * sqd_short,
                                      int seed, int global_len,
                                      sgeneration_mode_t mode);


/** 
    Makes one sequences longer given a model. See sgenerate_extensions for details.
    @return pointer to the whole sequence
    @param smo:        given model
    @param O:          given sequence to make longer
    @param len:        original length of sequence
    @param new_len:    wanted length of sequence
    @param alpha:
    @param mode:
*/
  double *sgenerate_single_ext (smodel * smo, double *O, const int len,
                                int *new_len, double **alpha,
                                sgeneration_mode_t mode);


/** Generate a single next value based on a trained model and on a seq of
   length "len". Use the most prob. state given the seq as an initial state
   and determin the next state und the symbol with the RNG. 
   @param smo:        given model
   @param O:          given sequence 
   @param len:        length of sequence
*/
  double sgenerate_next_value (smodel * smo, double *O, const int len);

/*@} sgenerate section */

#ifdef __cplusplus
}
#endif
#endif                          /* SGENERATE_H */
