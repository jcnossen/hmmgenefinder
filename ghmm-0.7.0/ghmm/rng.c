/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/rng.c
*       Authors:  Alexander Schliep
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
*       This file is version $Revision: 1276 $
*                       from $Date: 2005-08-24 17:34:29 +0200 (Wed, 24 Aug 2005) $
*             last change by $Author: cic99 $.
*
*******************************************************************************/

#ifdef WIN32
#  include "win_config.h"
#endif

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#include <stdlib.h>
#include <stdio.h>
#include "rng.h"
#include "math.h"
#include "time.h"

/* The global RNG */
GHMM_RNG *RNG;

#ifdef _MSC_VER
	#define SRANDOM(s) srand(s)
	#define RANDOM() rand()
#else
	#define SRANDOM(s) srandom(s)
	#define RANDOM() random()
#endif

#ifndef DO_WITH_GSL

ghmm_rng_state rng_state;
char rng_name[] = "random";

void ghmm_rng_set (GHMM_RNG * r, unsigned long int seed)
{
  SRANDOM(seed);
}

double ghmm_rng_uniform (GHMM_RNG * r)
{
  return ((double) RANDOM ()) / (RAND_MAX + 1.0);
}

const char *ghmm_rng_name (GHMM_RNG * r)
{
  return rng_name;
}

void ghmm_rng_init (void)
{
#ifndef _MSC_VER
  initstate (1, rng_state, sizeof (ghmm_rng_state));
#endif
  RNG = &rng_state;
}

#else

void ghmm_rng_init (void)
{
  gsl_rng_env_setup ();
  RNG = gsl_rng_alloc (gsl_rng_default);
}

#endif

void ghmm_rng_timeseed (GHMM_RNG * r)
{
  unsigned long tm;             /* Time seed */
  unsigned int timeseed;

  timeseed = time (NULL);
  srand (timeseed);
  tm = rand ();
  GHMM_RNG_SET (r, tm);
  /*printf("# using rng '%s' seed=%ld\n", GHMM_RNG_NAME(r), tm);  */
  fflush (stdout);
}



/* End of: rng.c */
