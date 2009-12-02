/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/sdfoba.h
*       Authors:  Wasinee Rungsarityotin, Utz Pape, Andrea Weisse
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
*       This file is version $Revision: 1264 $ 
*                       from $Date: 2005-08-10 18:33:54 +0200 (Wed, 10 Aug 2005) $
*             last change by $Author: grunau $.
*
*******************************************************************************/
#ifndef SDFOBA_H
#define SDFOBA_H

#ifdef __cplusplus
extern "C" {
#endif

/**@name HMM-Modell */
/*@{ (Doc++-Group: model) */

int sdfoba_forward (sdmodel * mo, const int *O, int len, double **alpha,
                    double *scale, double *log_p);

int sdfoba_backward (sdmodel * mo, const int *O, int len, double **beta,
                     const double *scale);

  int sdfoba_logp (sdmodel * mo, const int *O, int len, double *log_p);


#ifdef __cplusplus
}
#endif
#endif
