/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/matrix.c
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

#include <math.h>
#include <float.h>
#include "matrix.h"
#include "vector.h"
#include "model.h"
#include "rng.h"
#include "randvar.h"
#include <ghmm/mes.h>
#include <ghmm/internal.h>

static void lrdecomp (int dim, double **a, double *p);
static void lyequalsb (double **a, double *b, double *p, int dim, double *y);
static void ltranspxequalsy (double **a, double *y, double *p, int dim,
                             double *x);

/*============================================================================*/

int matrix_d_read (scanner_t * s, double **matrix, int max_zeile,
                   int max_spalte)
{
#define CUR_PROC "matrix_d_read"
  int len = 0, zeile = 0;
  scanner_consume (s, '{');
  if (s->err)
    return (-1);
  while (!s->eof && !s->err && s->c - '}') {
    if (zeile >= max_zeile) {
      scanner_error (s, "too many rows in matrix");
      return (-1);
    }
    /* Memory allocation takes place in scanner_get_double_earray */
    matrix[zeile] = scanner_get_double_earray (s, &len);
    if (len != max_spalte) {
      scanner_error (s, "wrong number of elements in matrix");
      return (-1);
    }
    scanner_consume (s, ';');
    if (s->err) {
      scanner_error (s, "missing ';' or wrong number of columns");
      return (-1);
    }
    zeile++;
  }
  scanner_consume (s, '}');
  if (s->err)
    return (-1);
  if (zeile < max_zeile) {
    scanner_error (s, "rows missing in matrix");
    return (-1);
  }
  return (0);
#undef CUR_PROC
}                               /* matrix_d_read */

/*============================================================================*/

int matrix_i_read (scanner_t * s, int **matrix, int max_zeile, int max_spalte)
{
#define CUR_PROC "matrix_i_read"
  int len = 0, zeile = 0;
  scanner_consume (s, '{');
  if (s->err)
    return (-1);
  while (!s->eof && !s->err && s->c - '}') {
    if (zeile >= max_zeile) {
      scanner_error (s, "too many rows in matrix");
      return (-1);
    }
    /* Memory allocation takes place in scanner_get_int_array */
    matrix[zeile] = scanner_get_int_array (s, &len);
    if (len != max_spalte) {
      scanner_error (s, "wrong number of elements in matrix");
      return (-1);
    }
    scanner_consume (s, ';');
    if (s->err) {
      scanner_error (s, "missing ';' or wrong number of columns");
      return (-1);
    }
    zeile++;
  }
  scanner_consume (s, '}');
  if (s->err)
    return (-1);
  if (zeile < max_zeile) {
    scanner_error (s, "rows missing in matrix");
    return (-1);
  }
  return (0);
#undef CUR_PROC
}                               /* matrix_i_read */

/*============================================================================*/

/* allocation of matrices with fixed dimensions  */
double **stat_matrix_d_alloc (int n, int m)
{
#define CUR_PROC "stat_matrix_d_alloc"
  int i;
  double **A;
  double *tmp;
 
  /* if (!(A = mes_calloc (n * sizeof (*A) + n * m * sizeof (**A)))) { */


  if (!(A = mes_calloc (n * sizeof (double*) + n * m * sizeof (double)))) {
    mes_proc ();
    goto STOP;
  }
  
  tmp = (double *) (A + n);
  for (i = 0; i < n; i++) {
    A[i] = tmp;
    tmp += m;
  }
  return A;
STOP:
  stat_matrix_d_free (&A);
  return NULL;
#undef CUR_PROC
}


int stat_matrix_d_free (double ***matrix)
{
#define CUR_PROC "stat_matrix_d_free"
  mes_check_ptr (matrix, return (-1));
  if (!*matrix)
    return (0);
  free (*matrix);
  return 0;
#undef CUR_PROC
}

/*============================================================================*/

/* allocation of matrices with fixed dimensions  */
int **stat_matrix_i_alloc (int n, int m)
{
#define CUR_PROC "stat_matrix_i_alloc"
  int i;
  int **A;
  int *tmp;
  if (!(A = mes_calloc (n * sizeof(int*) + n * m * sizeof (int)))) {
    mes_proc ();
    goto STOP;
  }

  tmp = (int *) (A + n);
  for (i = 0; i < n; i++) {
    A[i] = tmp;
    tmp += m;
  }
  return A;
STOP:
  stat_matrix_i_free (&A);
  return NULL;
#undef CUR_PROC
}


int stat_matrix_i_free (int ***matrix)
{
#define CUR_PROC "stat_matrix_i_free"
  mes_check_ptr (matrix, return (-1));
  if (!*matrix)
    return (0);
  free (*matrix);
  return 0;
#undef CUR_PROC
}


/*============================================================================*/


double **matrix_d_alloc (int zeilen, int spalten)
{
#define CUR_PROC "matrix_d_alloc"
  double **matrix;
  int i;

  /*printf("*** matrix_d_alloc %d zeilen, %d spalten:\n",zeilen, spalten);*/

  ARRAY_CALLOC (matrix, zeilen);
  for (i = 0; i < zeilen; i++)
    ARRAY_CALLOC (matrix[i], spalten);
  return matrix;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  matrix_d_free (&matrix, zeilen);
  return NULL;
#undef CUR_PROC
}                               /* matrix_d_alloc */


double *** matrix3d_d_alloc(int i, int j, int k) {
#define CUR_PROC "matrix3d_d_alloc"
  double *** matrix;
  int a, b;
  
  /* printf("*** matrix_d_alloc %d zeilen, %d spalten:\n",zeilen, spalten); */
  
  ARRAY_CALLOC (matrix, i);
  for (a = 0; a < i; a++) {
    ARRAY_CALLOC (matrix[a], j);
    for (b=0; b<j; b++) {
      ARRAY_CALLOC (matrix[a][b], k);
    }
  }
  return matrix;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  matrix3d_d_free(&matrix, i, j);
  return NULL;
#undef CUR_PROC
} /* matrix3d_d_alloc */

/** gets a pointer on a 3d matrix and rows and cols **/
int matrix3d_d_free(double **** matrix, int i, int j) {
# define CUR_PROC "matrix3d_d_free"
  int a,b;
  mes_check_ptr(matrix, return(-1));
  if ( !*matrix) return(0);
  for (a = i - 1; a >=  0; a--) {
    for (b=j-1; b>=0; b--)
      m_free((*matrix)[a][b]);
    m_free((*matrix)[a]);
  }
  m_free(*matrix);
  return (0);
# undef CUR_PROC
} /* matrix3d_d_free */

int matrix_d_free(double ***matrix, long zeilen) {
# define CUR_PROC "matrix_d_free"
  long i;
  mes_check_ptr (matrix, return (-1));
  if (!*matrix)
    return (0);
  for (i = zeilen - 1; i >= 0; i--)
    m_free ((*matrix)[i]);
  m_free (*matrix);
  return (0);
# undef CUR_PROC
}                               /* matrix_d_free */



/*============================================================================*/



double **matrix_d_alloc_copy (int zeilen, int spalten, double **copymatrix)
{
#define CUR_PROC "matrix_d_alloc_copy"
  double **matrix;
  int i, j;
  ARRAY_CALLOC (matrix, zeilen);
  for (i = 0; i < zeilen; i++) {
    ARRAY_CALLOC (matrix[i], spalten);
    for (j = 0; j < spalten; j++)
      matrix[i][j] = copymatrix[i][j];
  }
  return matrix;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  matrix_d_free (&matrix, zeilen);
  return NULL;
#undef CUR_PROC
}                               /* matrix_d_alloc_copy */

/*============================================================================*/

int **matrix_i_alloc (int zeilen, int spalten)
{
#define CUR_PROC "matrix_i_alloc"
  int **matrix;
  int i;
  ARRAY_CALLOC (matrix, zeilen);
  for (i = 0; i < zeilen; i++)
    ARRAY_CALLOC (matrix[i], spalten);
  return matrix;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  matrix_i_free (&matrix, zeilen);
  return NULL;
#undef CUR_PROC
}                               /* matrix_i_alloc */

/*============================================================================*/

int matrix_i_free (int ***matrix, long zeilen)
{
# define CUR_PROC "matrix_i_free"
  long i;
  mes_check_ptr (matrix, return (-1));
  if (!*matrix)
    return (0);
  for (i = 0; i < zeilen; i++) {
    m_free ((*matrix)[i]);
  }
  m_free (*matrix);
  return (0);
# undef CUR_PROC
}                               /* matrix_i_free */

/*============================================================================*/
/*============================================================================*/

int*** matrix3d_i_alloc(int zeilen, int spalten, int hoehe) {
#define CUR_PROC "matrix_i_alloc"
  int ***matrix;
  int i, j;
  ARRAY_CALLOC (matrix, zeilen);
  for (i = 0; i < zeilen; i++) {
    ARRAY_CALLOC (matrix[i], spalten);
    for (j=0; j < spalten; j++)
      ARRAY_CALLOC (matrix[i][j], hoehe);
  }
  return matrix;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  matrix3d_i_free(&matrix, zeilen, spalten);
  return NULL;
#undef CUR_PROC
} /* matrix3d_i_alloc */

/*============================================================================*/

int matrix3d_i_free(int **** matrix, int zeilen, int spalten) {
# define CUR_PROC "matrix_i_free"
  long i, j;
  mes_check_ptr(matrix, return(-1));
  if ( !*matrix ) return(0);
  for (i = zeilen - 1; i >= 0; i--) {
    for (j=spalten - 1; j >= 0; j--) 
      m_free((*matrix)[i][j]);
    m_free((*matrix)[i]);
  }
  m_free(*matrix);
  return (0);
# undef CUR_PROC
} /* matrix3d_i_free */

/*============================================================================*/
void matrix_d_print(FILE *file, double **matrix, int zeilen, int spalten, 
		    char *tab, char *separator, char *ending) {
  int i;
  for (i = 0; i < zeilen; i++)
    vector_d_print (file, matrix[i], spalten, tab, separator, ending);
}                               /* matrix_d_print */

/*============================================================================*/

void matrix_d_print_prec (FILE * file, double **matrix, int zeilen,
                          int spalten, int width, int prec, char *tab,
                          char *separator, char *ending)
{
  int i;
  for (i = 0; i < zeilen; i++)
    vector_d_print_prec (file, matrix[i], spalten, width, prec,
                         tab, separator, ending);
}                               /* matrix_d_print_prec */

/*============================================================================*/

void matrix_i_print (FILE * file, int **matrix, int zeilen, int spalten,
                     char *tab, char *separator, char *ending)
{
  int i;
  for (i = 0; i < zeilen; i++)
    vector_i_print (file, matrix[i], spalten, tab, separator, ending);
}                               /* matrix_i_print */


/*============================================================================*/
int matrix_d_notzero_columns (double **matrix, int row, int max_col)
{
  int i, count = 0;
  for (i = 0; i < max_col; i++)
    if (matrix[row][i])
      count++;
  return count;
}                               /* matrix_d_notzero_columns */

/*============================================================================*/
int matrix_d_notzero_rows (double **matrix, int col, int max_row)
{
  int i, count = 0;
  for (i = 0; i < max_row; i++)
    if (matrix[i][col])
      count++;
  return count;
}                               /* matrix_d__notzero_rows */

/*============================================================================*/
/* Scales the row vectors of a matrix to have the sum 1 */
int matrix_d_normalize (double **matrix, int rows, int cols)
{
#define CUR_PROC "matrix_d_normalize"
  int i;
  for (i = 0; i < rows; i++)
    if (vector_normalize (matrix[i], cols) == -1)
      mes (MES_WIN, "WARNING: sum row[%d] == 0!\n", i);
  /* return (-1); */
  return 0;
#undef CUR_PROC
}                               /* matrix_d_normalize */

/*============================================================================*/
void matrix_d_random_values (double **matrix, int rows, int cols,
                             double min, double max)
{
  int i, j;
  double interval;
  if (max < min) {
    min = 0.0;
    max = 1.0;
  }
  interval = max - min;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      matrix[i][j] = min + GHMM_RNG_UNIFORM (RNG) * interval;
}                               /* matrix_d_random_values */

/*============================================================================*/
/* Fixed value for final state */
void matrix_d_random_const_values (double **matrix, int rows, int cols,
                                   double min, double max, double c)
{
  int i, j;
  double interval;
  if (rows < 1) {
    mes (MES_WIN, "WARNING: rows = %d not allowed\n", rows);
    return;
  }
  if (max < min) {
    min = 0.0;
    max = 1.0;
  }
  interval = max - min;
  for (i = 0; i < rows - 1; i++)
    for (j = 0; j < cols; j++)
      matrix[i][j] = min + GHMM_RNG_UNIFORM (RNG) * interval;
  for (j = 0; j < cols; j++)
    matrix[rows - 1][j] = c;
}                               /* matrix_d_random_const_values */


/*============================================================================*/
void matrix_d_const_values (double **matrix, int rows, int cols, double c)
{
  int i, j;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      matrix[i][j] = c;
}                               /* matrix_d_const_values */

/*============================================================================*/
void matrix_d_random_left_right (double **matrix, int rows, int cols)
{
  int i, j;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      if (j == i || j == i + 1)
        matrix[i][j] = GHMM_RNG_UNIFORM (RNG);
      else
        matrix[i][j] = 0.0;
}                               /* matrix_d_random_values */

/*============================================================================*/
void matrix_d_left_right_strict (double **matrix, int rows, int cols)
{
  int i, j;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      if (j == i + 1)
        matrix[i][j] = 1.0;
      else
        matrix[i][j] = 0.0;
}                               /* matrix_d_left_right_strict */

/*============================================================================*/
int matrix_d_gaussrows_values (double **matrix, int rows, int cols,
                               double **mue, double u)
{
# define CUR_PROC "matrix_gaussrows_values"
  int res = -1;
  double *mean;
  int i, j;
  if (u <= 0.0) {
    mes_prot ("sigma^2 <= 0.0 not allowed\n");
    goto STOP;
  }
  if (*mue == NULL) {
    /* for each row, a random mean value mean[i] in (0, cols-1) */
    ARRAY_CALLOC (mean, rows);
    for (i = 0; i < rows; i++)
      mean[i] = GHMM_RNG_UNIFORM (RNG) * (cols - 1);
    /* for (i = 0; i < rows; i++) printf("%6.4f ", mean[i]); printf("\n"); */
    *mue = mean;
  }
  else
    /* Check, if the mean value is on the right interval */
    mean = *mue;
  for (i = 0; i < rows; i++) {
    /* Gauss-distribution around the mean value for each state. */
    for (j = 0; j < cols; j++) {
      matrix[i][j] = randvar_normal_density ((double) j, mean[i], u);
      if (matrix[i][j] == -1) {
        mes_proc ();
        goto STOP;
      }
      /* To avoid zero: (Cheap version!!) */
      if (matrix[i][j] < 0.0001)
        matrix[i][j] = 0.0001;  /* The Output has only 4 significant digits! */
    }
  }
  res = 0;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return (res);
# undef CUR_PROC
}                               /* matrix_gaussrows_values */

/*============================================================================*/
void matrix_d_const_preserve_struct (double **matrix, int rows, int cols,
                                     double c)
{
  int i, j;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++) {
      if (matrix[i][j] != 0)
        matrix[i][j] = c;
    }
}                               /* matrix_d_const_preserve_struct */

/*============================================================================*/
void matrix_d_random_preserve_struct (double **matrix, int rows, int cols)
{
  int i, j;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++) {
      if (matrix[i][j] != 0)
        matrix[i][j] = GHMM_RNG_UNIFORM (RNG);
    }
}                               /* matrix_d_random_preserve_struct */

/*============================================================================*/
void matrix_d_transpose (double **A, int rows, int cols, double **A_T)
{
  int i, j;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      A_T[j][i] = A[i][j];
}                               /* matrix_d_transpose */


/*----------------------------------------------------------------------------*/
/* Decomposes a positiv definit, symmetric matrix A in L*L-T, that is except for
   the diagonal, L (= lower triangular marix) is stored in A. The Diagonal is 
   stored in a vector p. 
*/
static void lrdecomp (int dim, double **a, double *p)
{
  int k, i, j;
  double x;
  for (i = 0; i < dim; i++) {
    for (j = i; j < dim; j++) {
      x = a[i][j];
      for (k = i - 1; k >= 0; k--)
        x = x - a[j][k] * a[i][k];
      if (i == j) {
        if (x < DBL_MIN)
          mes (MES_WIN, "FEHLER: Pivotel.<=0!");
        p[i] = 1 / sqrt (x);
      }
      else
        a[j][i] = x * p[i];
    }
  }
}

/*----------------------------------------------------------------------------*/
/* Solves L*y=b, L is a lower triangular matrix stored in A, p=1/diagonal elements. */
static void lyequalsb (double **a, double *b, double *p, int dim, double *y)
{
  int k, j;
  for (k = 0; k < dim; k++) {
    y[k] = b[k];
    for (j = 0; j < k; j++)
      y[k] = y[k] - a[k][j] * y[j];
    y[k] = y[k] * p[k];
  }
}

/*----------------------------------------------------------------------------*/
/* Solves L-T*x=y, L-T an upper triangular matrix, BUT saved in A as L! */
static void ltranspxequalsy (double **a, double *y, double *p, int dim,
                             double *x)
{
  int k, j;
  for (k = dim - 1; k >= 0; k--) {
    x[k] = y[k];
    for (j = k + 1; j < dim; j++)
      x[k] = x[k] - a[j][k] * x[j];
    x[k] = x[k] * p[k];
  }
}


/*============================================================================*/
/* Solves a linear equation system for a symmetric, positiv definite matrix. */
int matrix_cholesky (double **a, double *b, int dim, double *x)
{
#define CUR_PROC "matrix_cholesky"
  int res = -1;
  double *p, *y;
  ARRAY_CALLOC (p, dim);
  ARRAY_CALLOC (y, dim);
  lrdecomp (dim, a, p);
  lyequalsb (a, b, p, dim, y);
  ltranspxequalsy (a, y, p, dim, x);
  res = 0;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return (res);
#undef CUR_PROC
}

/* Finds the determinant of a symetric, positiv definit matrix. */
int matrix_det_symposdef (double **a, int dim, double *det)
{
#define CUR_PROC "matrix_det_symposdef"
  int res = -1;
  int i;
  double *p, r;
  ARRAY_CALLOC (p, dim);
  lrdecomp (dim, a, p);
  *det = 1.0;
  for (i = 0; i < dim; i++) {
    r = 1 / p[i];
    *det *= (r * r);
  }
  res = 0;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return (res);
#undef CUR_PROC
}

/*============================================================================*/
/* Copies a matrix, the memory allocations must to be done outside! */
void matrix_d_copy (double **src, double **target, int rows, int cols)
{
  int i, j;

  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      target[i][j] = src[i][j];
}

/*============================================================================*/
/**
  Checks whether a quadratic double matrix is stochastic
  @return 0/1 flag for true/false
  @param  double NxN matrix to be checked
  @param  matrix dimension N (matrix must be quadaratic)
  */
int matrix_d_check_stochasticity (double **matrix, int N)
{
  int i, j;
  double row_sum;
  int stochastic = 1;

  for (i = 0; i < N; i++) {
    row_sum = 0.0;
    for (j = 0; j < N; j++) {
      row_sum = row_sum + matrix[i][j];
    }
    if (row_sum != 1.0) {
      stochastic = 0;
      break;
    }
  }
  return (stochastic);
}
