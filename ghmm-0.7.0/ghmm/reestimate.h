/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/reestimate.h
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
*       This file is version $Revision: 1355 $ 
*                       from $Date: 2005-09-14 13:57:38 +0200 (Wed, 14 Sep 2005) $
*             last change by $Author: cic99 $.
*
*******************************************************************************/
#ifndef REESTIMATE_H
#define REESTIMATE_H
#include <ghmm/sequence.h>
#include <ghmm/model.h>

#ifdef __cplusplus
extern "C" {
#endif


/** matrix allocation and free for training algorithms */
  int reestimate_alloc_matvek (double ***alpha, double ***beta,
                               double **scale, int T, int N);
  int reestimate_free_matvek (double **alpha, double **beta, double *scale,
                              int T);


/**@name Baum-Welch-Algorithm */
/*@{ (Doc++-Group: reestimate) */

/** Baum-Welch-Algorithm for parameter reestimation (training) in
    a discrete (discrete output functions) HMM. Scaled version
    for multiple sequences, alpha and beta matrices are allocated with
	stat_matrix_d_alloc 
    New parameters set directly in hmm (no storage of previous values!).
    For reference see:  
    Rabiner, L.R.: "`A Tutorial on Hidden {Markov} Models and Selected
                Applications in Speech Recognition"', Proceedings of the IEEE,
	77, no 2, 1989, pp 257--285    
  @return            0/-1 success/error
  @param mo          initial model
  @param sq          training sequences
  */

  int reestimate_baum_welch (model * mo, sequence_t * sq);

/** Just like reestimate_baum_welch, but you can limit
    the maximum number of steps
  @return            0/-1 success/error
  @param mo          initial model
  @param sq          training sequences
  @param max_step    maximal number of Baum-Welch steps
  @param likelihood_delta minimal improvement in likelihood required for carrying on. Relative value
  to log likelihood
  */
  int reestimate_baum_welch_nstep (model * mo, sequence_t * sq, int max_step,
                                   double likelihood_delta);


/** Update the emissions according to the tie groups by computing the mean
    values within all groups.
    */
  void reestimate_update_tie_groups (model * mo);


/** Baum-Welch-Algorithm for parameter reestimation (training) in
    a StateLabelHMM. Scaled version for multiple sequences, alpha and 
    beta matrices are allocated with stat_matrix_d_alloc 
    New parameters set directly in hmm (no storage of previous values!).
    For reference see:  
    Rabiner, L.R.: "`A Tutorial on Hidden {Markov} Models and Selected
                Applications in Speech Recognition"', Proceedings of the IEEE,
	77, no 2, 1989, pp 257--285    
  @return            0/-1 success/error
  @param mo          initial model
  @param sq          training sequences
  */

  int reestimate_baum_welch_label (model * mo, sequence_t * sq);

/** Just like reestimate_baum_welch_label, but you can limit
    the maximum number of steps
  @return                   0/-1 success/error
  @param mo                 initial model
  @param sq                 training sequences
  @param max_step           maximal number of Baum-Welch steps
  @param likelihood_delta   minimal improvement in likelihood required for
                            carrying on. Relative value to log likelihood
  */
  int reestimate_baum_welch_nstep_label (model * mo, sequence_t * sq,
                                         int max_step,
                                         double likelihood_delta);



#ifdef __cplusplus
}
#endif
#endif                          /* REESTIMATE_H */
/*@} (Doc++-Group: reestimate) */
