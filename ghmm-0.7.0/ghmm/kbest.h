/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/kbest.h
*       Authors:  Alexander Riemer, Janne Grunau
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
#ifndef KBEST_H
#define KBEST_H

#ifdef __cplusplus
extern "C" {
#endif

/**
   Calculates the most probable labeling for the given sequence in the given
   model using k-best decoding.
   Labels must be from interval [0:max_label] without gaps!!! (not checked)
   Model must not have silent states. (checked in Python wrapper)
   @return array of labels (internal representation)
   @param mo:         pointer to a model
   @param o_seq:      output sequence (array of internal representation chars)
   @param seq_len:    length of output sequence
   @param k:          number of hypotheses to keep for each state
   @param log_p:      variable reference to store the log prob. of the labeling
 */
  int *kbest (model * mo, int *o_seq, int seq_len, int k, double *log_p);

#ifdef __cplusplus
}
#endif
#endif
