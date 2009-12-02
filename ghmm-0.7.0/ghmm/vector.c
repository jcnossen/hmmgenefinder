/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/vector.c
*       Authors:  Bernhard Knab
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

#ifdef WIN32
#  include "win_config.h"
#endif

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <float.h>
#include <stdio.h>
#include "mes.h"
#include "vector.h"
#include "rng.h"


/*============================================================================*/
/* Scales the elements of a vector to have the sum 1 */
/* PROBLEM: Entries can get very small and be rounded to 0 */
int vector_normalize (double *v, int len)
{
#define CUR_PROC "vector_normalize"
  int i;
  double sum = 0.0;
  for (i = 0; i < len; i++)
    sum += v[i];
  if (sum < DBL_MIN) {
    mes_prot ("Can't normalize vector. Sum eq. zero \n");
    return (-1);
  }
  for (i = 0; i < len; i++)
    v[i] /= sum;
  return 0;
#undef CUR_PROC
}                               /* vector_normalize */


/*============================================================================*/
void vector_const_values (double *v, int len, double c)
{
  int i;
  for (i = 0; i < len; i++)
    v[i] = c;
}                               /* vector_const_values */


/*============================================================================*/
void vector_const_preserve_struct (double *v, int len, double c)
{
  int i;
  for (i = 0; i < len; i++) {
    if (v[i] != 0.0)
      v[i] = c;
  }
}                               /* vector_const_preserve_struct */


/*============================================================================*/
void vector_random_values (double *v, int len)
{
  int i;
  for (i = 0; i < len; i++)
    v[i] = GHMM_RNG_UNIFORM (RNG);
}                               /* vector_random_values */


/*============================================================================*/
void vector_random_preserve_struct (double *v, int len)
{
  int i;
  for (i = 0; i < len; i++) {
    if (v[i] != 0.0)
      v[i] = GHMM_RNG_UNIFORM (RNG);
  }
}                               /* vector_random_preserve_struct */


/*============================================================================*/
void vector_d_print (FILE * file, double *vector, int len,
                     char *tab, char *separator, char *ending)
{
  int j;

  fprintf (file, "%s", tab);
  if (len > 0)
    fprintf (file, "%6.2f", vector[0]);

  for (j = 1; j < len; j++)
    fprintf (file, "%s %6.2f", separator, vector[j]);
  fprintf (file, "%s\n", ending);
}                               /* vector_d_print */

/*============================================================================*/
void vector_d_print_prec (FILE * file, double *vector, int len, int width,
                          int prec, char *tab, char *separator, char *ending)
{
  int j;
  if (len > 0)
    fprintf (file, "%s%*.*f", tab, width, prec, vector[0]);
  for (j = 1; j < len; j++)
    fprintf (file, "%s %*.*f", separator, width, prec, vector[j]);
  fprintf (file, "%s\n", ending);
}                               /* vector_d_print */

/*============================================================================*/
void vector_i_print (FILE * file, int *vector, int len,
                     char *tab, char *separator, char *ending)
{
  int j;
  fprintf (file, "%s", tab);
  if (len > 0)
    fprintf (file, "%3d", vector[0]);
  for (j = 1; j < len; j++)
    fprintf (file, "%s %3d", separator, vector[j]);
  fprintf (file, "%s\n", ending);
}                               /* vector_i_print */

/*============================================================================*/
int vector_mat_times_vec (double **A, double *x, int n, int m, double *v)
{
#define CUR_PROC "vector_mat_times_vec"
  int i, j;
  for (i = 0; i < n; i++) {
    v[i] = 0.0;
    for (j = 0; j < m; j++)
      v[i] += A[i][j] * x[j];
  }
  return 0;
#undef CUR_PROC
}                               /* vector_mat_times_vec */
