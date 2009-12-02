/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/const.h
*       Authors:  Bernhard Knab, Benjamin Georgi
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
*       This file is version $Revision: 1364 $
*                       from $Date: 2005-09-16 12:53:05 +0200 (Fri, 16 Sep 2005) $
*             last change by $Author: cic99 $.
*
*******************************************************************************/

#ifndef CONST_H
#define CONST_H

/**@name GHMM Constants */
/*@{ (Doc++-Group: Constants) */


/** 
    Convergence: Halt criterium for Baum-Welch reestimation if the difference
    of log(P) in two consecutive iterations is smaller than (EPS\_ITER\_BW * log(P))..
*/
/* #define EPS_ITER_BW      0.00001 */
#define EPS_ITER_BW      0.0001

/**
  If the absolute difference of two numbers is smaller the EPS_PREC, then the
  numbers are equal. (Instead of using zero )
  */
#define EPS_PREC         1E-8

/**
  Minimum value for U
  */
#define EPS_U            1E-4

/**
  Maximum value for U (Turning point at 100 ?)
  */
#define MAX_U            10000

/**
  Maximum number of iterations in reestimate
  */
#define MAX_ITER_BW      500

/**
  Maximum length of a sequence
  */
#define MAX_SEQ_LEN 100000
/* #define MAX_SEQ_LEN      1000000 */

/**
  Maximum number of sequences 
  */
#define MAX_SEQ_NUMBER   1500000

/* in float.h: DBL_EPSILON = 0.000000000000000222044604925031... 
*/

#if 0
/* Now a member of smodel */
#define COS 1
#endif

/**
   A value that is put in for log_p in the calculation of
   the objective function if sfoba_logp returns -1 (error).
*/
#define PENALTY_LOGP -500.0

/**
   The left limit for the normal density
*/
#define EPS_NDT  0.1

/*@} (Doc++-Group: Constants) */

#endif /* CONST_H */
