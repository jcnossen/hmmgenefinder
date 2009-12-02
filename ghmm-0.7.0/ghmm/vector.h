/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/vector.h
*       Authors:  Bernhard Knab, Peter Pipenbacher
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

#ifndef VECTOR_H
#define VECTOR_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

/**@name Vektor-Funktionen */
/*@{ (Doc++-Group: vector) */

/**
  Scales the sum of elements in a vector to one.
  @return 0 for success; -1 for error
  @param v    vector
  @param len  length of the vector       
*/
  int vector_normalize (double *v, int len);

/**
  Gives all elements in a vector a constant value
  @param v    vector
  @param len  length of the vector
  @param c    given value for the elements
*/
  void vector_const_values (double *v, int len, double c);

/**
  Gives all elements, not equal zero, in a vector a constant value
  @param v    vector
  @param len  length of the vector
  @param c    given value for the elements
*/
  void vector_const_preserve_struct (double *v, int len, double c);

/**
  Gives all elements in a vector random values between 0 and 1
  @param v    vector
  @param len  length of the vector       
*/
  void vector_random_values (double *v, int len);

/**
  Gives all elements, not equal zero, in a vector random values between 0 and 1
  @param v    vector
  @param len  length of the vector   
*/
  void vector_random_preserve_struct (double *v, int len);

/**
  Writes a double vector (without parenthesis)
  @param file       output file
  @param vector     vector to write
  @param len        dimension
  @param tab        format: leading tabs
  @param separator  format: separator for columns
  @param ending     format: end of a row  
  */
  void vector_d_print (FILE * file, double *vector, int len,
                       char *tab, char *separator, char *ending);

/**
  Writes a double vector (without parenthesis) with given number of decimal places
  @param file       output file
  @param vector     vector to write
  @param len        dimension
  @param width      format: total number of decimal places
  @param prec       format: number of decimal places after the comma
  @param tab        format: leading tabs
  @param separator  format: separator for columns
  @param ending     format: end of a row 
  */
  void vector_d_print_prec (FILE * file, double *vector, int len, int width,
                            int prec, char *tab, char *separator,
                            char *ending);

/**
  Writes an integer vector (without parenthesis)
  @param file       output file
  @param vector     vector to write
  @param len        dimension
  @param tab        format: leading tabs
  @param separator  format: separator for columns
  @param ending     format: end of a row  
  */
  void vector_i_print (FILE * file, int *vector, int len,
                       char *tab, char *separator, char *ending);
/**
  Calculates Ax, where A is a double matrix and x a double vector
  @param A       n x m matrix
  @param x       vector to calculate
  @param n       number of rows
  @param m       number of columns
  @param v       calculated vector (return value)
  */
  int vector_mat_times_vec (double **A, double *x, int n, int m, double *v);

#ifdef __cplusplus
}
#endif
#endif                          /* VECTOR_H */
/*@} (Doc++-Group: vector) */
