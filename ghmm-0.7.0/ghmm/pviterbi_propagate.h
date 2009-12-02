/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/pviterbi_propagate.h
*       Authors:  Matthias Heinig, Janne Grunau
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
*       This file is version $Revision: 1328 $
*                       from $Date: 2005-09-09 17:18:54 +0200 (Fri, 09 Sep 2005) $
*             last change by $Author: grunau $.
*
*******************************************************************************/

#ifndef PVITERBI_PROPAGATE_H
#define PVITERBI_PROPAGATE_H
#ifdef __cplusplus
extern "C" {
#endif

#include "pmodel.h"
#include "psequence.h"
#include "linkedlist.h"

/*------------        Here comes the Propagate stuff          ------------- */

int * pviterbi_propagate(pmodel *mo, psequence * X, psequence * Y,
			 double *log_p, int *path_length, double max_size);


int * pviterbi_propagate_segment (pmodel *mo, psequence * X, psequence * Y,
				  double *log_p, int *path_length,
				  double max_size, int start_x, int start_y,
				  int stop_x, int stop_y, int start_state,
				  int stop_state, double start_log_p,
				  double stop_log_p);

#ifdef __cplusplus
}
#endif

#endif
