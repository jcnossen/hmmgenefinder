/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/smixturehmm.h
*       Authors:  Bernd Wichern
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
#ifndef SMIXTUREHMM_H
#define SMIXTUREHMM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <ghmm/sequence.h>
#include <ghmm/smodel.h>

/**
   @name mixture-shmm methods
 */

/*@{ smixturehmm section
 */

/**
 */
  int smixturehmm_cluster (FILE * outfile, double **cp, sequence_d_t * sqd,
                           smodel ** smo, int smo_number);

/**
 */
  double smixturehmm_like (smodel ** smo, int smo_number,
                           sequence_d_t * sqd_test, long *errors);

/**
 */
  int smixturehmm_init (double **cp, sequence_d_t * sqd, smodel ** smo,
                        int smo_number, int mode);

/**
 */
  int smixturehmm_calc_priors (double **cp, sequence_d_t * sqd, smodel ** smo,
                               int smo_number);

/**
 */
  int smixturehmm_calc_cp (double **cp, sequence_d_t * sqd, smodel ** smo,
                           int smo_number, double *total_train_w);

/**
*/
  void smixture_calc_logp (double **logp, int **error, sequence_d_t * sqd,
                           smodel ** smo, int smo_number);

/**
 */
  void smixturehmm_print_header (FILE * file, char *argv[], int flag);


/**
 */
  double *smixturehmm_avg_like (double **cp, sequence_d_t * sqd,
                                smodel ** smo, int smo_number);
#ifdef __cplusplus
}
#endif
/*@} smixturehmm section */
#endif                          /* SMIXTUREHMM_H */
