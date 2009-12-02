/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/gradescent.h
*       Authors:  Janne Grunau, Alexander Riemer, Benjamin Georgi
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
*       This file is version $Revision: 1294 $ 
*                       from $Date: 2005-09-02 19:32:52 +0200 (Fri, 02 Sep 2005) $
*             last change by $Author: grunau $.
*
*******************************************************************************/
#ifndef GRADESCENT_H
#define GRADESCENT_H

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------------*/
/**
   Trains the model with a set of annotated sequences till convergence using
   gradient descent.
   Model must not have silent states. (checked in Python wrapper)
   @return            0/-1 success/error
   @param mo:         pointer to a model
   @param sq:         struct of annotated sequences
   @param eta:        intial parameter eta (learning rate)
   @param no_steps    number of training steps
 */
int gradient_descent (model ** mo, sequence_t * sq, double eta, int no_steps);

/*----------------------------------------------------------------------------*/
/**
   computes matrices of n and m variables (expected values for how often a
   certain parameter from A or B is used)
   computes Baum-Welch variables implicit 
   @return                 0/-1 success/error
   @param mo:              pointer to a model
   @param alpha:           matrix of forward variables
   @param backward:        matrix of backward variables
   @param scale:           scaling vector from forward-backward-algorithm
   @param seq:             sequence in internal representation
   @param seq_len:         length of sequence
   @param matrix_b:        matrix for parameters from B (n_b or m_b)
   @param matrix_a:        matrix for parameters from A (n_a or m_a)
   @param vec_pi:          vector for parameters in PI (n_pi or m_pi)
 */
int gradescent_compute_expectations (model * mo, double **alpha, double **beta,
		                     double *scale, int *seq, int seq_len,
				     double **matrix_b, double *matrix_a,
				     double *vec_pi);


#ifdef __cplusplus
}
#endif
#endif
