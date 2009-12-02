/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/root_finder.h
*       Authors:  Achim Gaedke
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
#ifndef ROOT_FINDER_H
#define ROOT_FINDER_H

#ifdef __cplusplus
extern "C" {
#endif
/**
   @name root finder
 */

/*@{ */

/**
   brent root finding algorithm.
   wrapps this functioncall to the gsl 1D-root finding interface
   @author Achim G\"adke
   @param x1 lower bracket value
   @param x2 upper bracket value
   @param tol tolerance for iteration
   @param A 1st extra parameter
   @param B 2nd extra parameter
   @param eps 3rd extra parameter
   @return root
 */
  double zbrent_AB (double (*func) (double, double, double, double),
                    double x1, double x2, double tol, double A, double B,
                    double eps);

/*@} */
#ifdef __cplusplus
}
#endif
#endif
