/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/mprintf.h
*       Authors:  Frank Nuebel
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
/*******************************************************************************
  author       : Frank Nübel
  filename     : ghmm/ghmm/mprintf.h
  created      : TIME: 11:27:32     DATE: Wed 14. May 1997
  $Id: mprintf.h 1191 2005-06-21 09:56:12Z cic99 $

Copyright (C) 1998-2005 Alexander Schliep
Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
Copyright (C) 2002-2005 Max-Planck-Institut fuer Molekulare Genetik, Berlin

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA




*******************************************************************************/

#ifndef MPRINTF_H
#define MPRINTF_H

#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif /*__cplusplus*/

/**
   @name Help functions for printing.
*/

  /**
   */
  char *mprintf (char *dst, int maxlen, char *format, ...);
  /**
   */
  char *mprintf_dyn (char *dst, int maxlen, char *format, ...);
  /**
   */
  char *mprintf_va (char *dst, int maxlen, char *format, va_list args);
  /**
   */
  char *mprintf_va_dyn (char *dst, int maxlen, char *format, va_list args);



#ifdef __cplusplus
}
#endif /*__cplusplus*/
#endif                          /* MPRINTF_H */
