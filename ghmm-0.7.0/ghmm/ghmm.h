/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/ghmm.h
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
*       This file is version $Revision: 1191 $
*                       from $Date: 2005-06-21 11:56:12 +0200 (Tue, 21 Jun 2005) $
*             last change by $Author: cic99 $.
*
*******************************************************************************/

#ifndef GHMM_H
#define GHMM_H

#ifdef __cplusplus
extern "C" {
#endif
/**@name GHMM-Globals */
/*@{ (Doc++-Group: globals) */
#ifdef __cplusplus
}
#endif
/** @name type_constants
    Constants giving model variations *//** Model is a left-right */
#define kNotSpecified (0)
#define kLeftRight (1)
/** Model contains silent states (i.e., states without emissions) */
#define kSilentStates (1 << 2)
/** Model has states with tied emission probabilities */
#define kTiedEmissions (1 << 3)
#define kUntied -1
/** Model has states emission probabilities conditioned on previous orders */
#define kHigherOrderEmissions (1 << 4)
#define kHasBackgroundDistributions (1 << 5)
#define kNoBackgroundDistribution -1
/** Model is a class HMM with labeled states */
#define kLabeledStates (1 << 6)
#endif
/*@} (Doc++-Group: GHMM-Globals) */
