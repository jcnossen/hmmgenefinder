/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/rng.h
*       Authors:  Alexander Schliep, Ben Rich
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
*       This file is version $Revision: 1268 $ 
*                       from $Date: 2005-08-17 23:24:06 +0200 (Wed, 17 Aug 2005) $
*             last change by $Author: grunau $.
*
*******************************************************************************/
#ifndef RNG_H
#define RNG_H

#ifndef DO_WITH_GSL

#include <stdlib.h>

typedef char ghmm_rng_state[8];
typedef ghmm_rng_state GHMM_RNG;

/* Functions */
#define GHMM_RNG_SET ghmm_rng_set
#define GHMM_RNG_UNIFORM ghmm_rng_uniform
#define GHMM_RNG_NAME ghmm_rng_name

#else

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Types defined in GSL */
typedef gsl_rng GHMM_RNG;

/* Functions defined in GSL */
#define GHMM_RNG_SET gsl_rng_set
#define GHMM_RNG_UNIFORM gsl_rng_uniform
#define GHMM_RNG_NAME gsl_rng_name

#endif /*  DO_WITH_GSL  */


/** @name rng Initialization for the random number generator */

/**
 */
extern GHMM_RNG *RNG;


#ifdef __cplusplus
extern "C" {
#endif

/**
 */
  void ghmm_rng_init (void);

/**
 */
  void ghmm_rng_timeseed (GHMM_RNG * r);


#ifndef DO_WITH_GSL
  void ghmm_rng_set (GHMM_RNG * aState, unsigned long int seed);
  double ghmm_rng_uniform (GHMM_RNG * r);
  const char *ghmm_rng_name (GHMM_RNG * r);
#endif

#ifdef __cplusplus
}
#endif
#endif                          /* RNG_H */
