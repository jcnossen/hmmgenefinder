/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/sreestimate.h
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
*       This file is version $Revision: 1206 $ 
*                       from $Date: 2005-07-05 22:19:51 +0200 (Tue, 05 Jul 2005) $
*             last change by $Author: grunau $.
*
*******************************************************************************/
#ifndef SREESTIMATE_H
#define SREESTIMATE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <ghmm/smodel.h>

/**@name SHMM-Baum-Welch-Algorithm */
/*@{ (Doc++-Group: sreestimate) */

/** Baum-Welch-Algorithm for parameter reestimation (training) in
    a continuous (continuous output functions) HMM. Scaled version
    for multiple sequences. Sequences may carry different weights 
    For reference see:  
    Rabiner, L.R.: "`A Tutorial on Hidden {Markov} Models and Selected
                Applications in Speech Recognition"', Proceedings of the IEEE,
	77, no 2, 1989, pp 257--285    
*/

/** @name struct smosqd_t
    structure that combines a continuous model (smo) and an integer
    sequence struct. Is used by sreestimate\_baum\_welch for 
    parameter reestimation.
 */
  struct smosqd_t {
  /** pointer of continuous model*/
    smodel *smo;
  /** sequence\_d\__t pointer */
    sequence_d_t *sqd;
  /** calculated log likelihood */
    double *logp;
  /** leave reestimation loop if diff. between successive logp values 
      is smaller than eps */
    double eps;
  /** max. no of iterations */
    int max_iter;
  };
  typedef struct smosqd_t smosqd_t;

  typedef struct local_store_t {
    int cos;
    double *pi_num;
    double pi_denom;
    double ***a_num;
    double **a_denom;
    double **c_num;
    double *c_denom;
    double **mue_num;
    double **u_num;
    double **mue_u_denom;       /* mue-denom. = u-denom. for sym. normal density */
    double **sum_gt_otot;       /* for truncated normal density */
    double **sum_gt_logb;       /* Control according to Q-function */
  } local_store_t;


/**
  Baum-Welch Algorithm for SHMMs.
  Training of model parameter with multiple double sequences (incl. scaling).
  New parameters set directly in hmm (no storage of previous values!). Matrices
  are allocated with stat_matrix_d_alloc.
  @return            0/-1 success/error
  @param cs         initial model and train sequences
  */
  int sreestimate_baum_welch (smosqd_t * cs);


#ifdef __cplusplus
}
#endif
#endif                          /* SREESTIMATE_H */
/*@} (Doc++-Group: sreestimate) */
