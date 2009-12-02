/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/foba.h
*       Authors:  Bernd Wichern, Benjamin Georgi
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
*       This file is version $Revision: 1323 $
*                       from $Date: 2005-09-09 13:31:33 +0200 (Fri, 09 Sep 2005) $
*             last change by $Author: grunau $.
*
*******************************************************************************/


#ifndef FOBA_H
#define FOBA_H

#include <ghmm/model.h>


#ifdef __cplusplus
extern "C" {
#endif

/**@name Forward-Backward-Algorithm 
   
   Forward-Backward-Algorithm for multiple integer
   sequences with scaling.
   For reference see:  
   Rabiner, L.R.: "`A Tutorial on Hidden {Markov} Models and Selected
                Applications in Speech Recognition"', Proceedings of the IEEE,
	77, no 2, 1989, pp 257--285
       
*/

/*@{ (Doc++-Group: foba) */

/** Forward-Algorithm.
  Calculates alpha[t][i], scaling factors scale[t] and log( P(O|lambda) ) for
  a given double sequence and a given model.
  @param mo:      model
  @param O:       sequence
  @param len:     length of sequence
  @param alpha:   alpha[t][i]
  @param scale:   scale factors
  @param log\_p:  log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
  int foba_forward (model * mo, const int *O, int length, double **alpha,
                    double *scale, double *log_p);

/** 
  Backward-Algorithm. 
  Calculates beta[t][i] given an integer sequence and a model. Scale factors 
  given as parameter (come from foba\_forward).
  @param mo:      model
  @param O:       sequence
  @param len:     length of sequence
  @param beta:    empty beta matrix
  @param scale:   scale factors
  @return 0 for success, -1 for error
  */
  int foba_backward (model * mo, const int *O, int len, double **beta,
                     const double *scale);

/** 
  Termination of Backward-Algorithm. 
  Calculates Backward-probability given an integer sequence, a model and
  the beta matrix. Scale factors given as parameter (come from foba\_forward).
  @param mo:      pointer to the model
  @param O:       sequence
  @param len:     length of sequence
  @param beta:    beta matrix
  @param scale    scale factors
  @param log\_p:  log probability
  @return 0 for success, -1 for error
  */
  int foba_backward_termination (model *mo, const int *O, int len,
				 double **beta, double *scale, double *log_p);

/**
  Calculation of  log( P(O|lambda) ). 
  Done by calling sfoba\_forward. Use this function if only the
  log likelihood and not alpha[t][i] is needed, alpha is allocated with
  stat_matrix_d_alloc
  @param  mo      model
  @param O        sequence
  @param len       length of sequence
  @param log\_p    log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
  int foba_logp (model * mo, const int *O, int len, double *log_p);


/** Forward-Algorithm (lean version).
  Calculates log( P(O|lambda) ) for a given double sequence and a given model.
  @param smo      model
  @param O        sequence
  @param length: length of sequence
  @param log\_p:  log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
  int foba_forward_lean (model * mo, const int *O, int len, double *log_p);


/* Labeled HMMs */

  int foba_label_forward (model * mo, const int *O, const int *label, int len,
                          double **alpha, double *scale, double *log_p);
  int foba_label_logp (model * mo, const int *O, const int *label, int len,
                       double *log_p);


  int foba_label_backward (model * mo, const int *O, const int *label,
                           int len, double **beta, double *scale,
                           double *log_p);



  int foba_initforward (model * mo, double *alpha_1, int symb, double *scale);
  double foba_stepforward (state * s, double *alpha_t, const double b_symb);

/*@} (Doc++-Group: foba) */

#ifdef __cplusplus
}
#endif
#endif
