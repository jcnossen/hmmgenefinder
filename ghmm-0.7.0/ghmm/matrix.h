/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/matrix.h
*       Authors:  Bernhard Knab, Benjamin Georgi
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
*       This file is version $Revision: 1232 $ 
*                       from $Date: 2005-07-27 18:34:17 +0200 (Wed, 27 Jul 2005) $
*             last change by $Author: ties $.
*
*******************************************************************************/
#ifndef MATRICE_H
#define MATRICE_H

#include <stdio.h>
#include <ghmm/scanner.h>

#ifdef __cplusplus
extern "C" {
#endif

/**@name Matrix */
/*@{ (Doc++-Group: matrix) */

/**
  Allocation of a double matrix. 
  @return pointer to a matrix
  @param rows: number of rows
  @param columns: number of columns
  */
  double **matrix_d_alloc (int n, int m);


/**
  Allocation of a 3 dimensional double matrix. 
  @return pointer to a matrix
  @param rows: number of rows
  @param columns: number of columns
  @param height: 3rd dimension
  */
double *** matrix3d_d_alloc(int i, int j, int k);
int matrix3d_d_free(double **** matrix, int i, int j);

int *** matrix3d_i_alloc(int i, int j, int k);
int matrix3d_i_free(int **** matrix, int i, int j);
/**
  Copying and allocation of a double matrix.
  @return pointer to a matrix
  @param rows: number of rows
  @param columns: number of columns
  @param copymatrix: matrix to copy 
  */
  double **matrix_d_alloc_copy (int rows, int columns, double **copymatrix);

/**
  Free the memory of a double matrix.
  @return 0 for succes; -1 for error
  @param  matrix: matrix to free
  @param  rows: number of rows
  */
  int matrix_d_free (double ***matrix, long zeilen);

/**
  Allocation of a static double matrix with a single malloc. 
  @return pointer to a matrix
  @param rows: number of rows
  @param columns: number of columns
  */
  double **stat_matrix_d_alloc (int n, int m);

/**
  Free the memory of a static double matrix.
  @return 0 for succes; -1 for error
  @param  matrix: matrix to free
  @param  rows: number of rows
  */
  int stat_matrix_d_free (double ***matrix);

/**
  Allocation of a static int matrix with a single malloc. 
  @return pointer to a matrix
  @param rows: number of rows
  @param columns: number of columns
  */
  int **stat_matrix_i_alloc (int n, int m);

/**
  Free the memory of a static int matrix.
  @return 0 for succes; -1 for error
  @param  matrix: matrix to free
  @param  rows: number of rows
  */
  int stat_matrix_i_free (int ***matrix);


/**
  Allocation of a integer matrix.
  @return pointer to a matrix
  @param rows: number of rows
  @param columns: number of columns
  */
  int **matrix_i_alloc (int rows, int columns);

/**
  Free the memory of a integer matrix.
  @return 0 for succes; -1 for error
  @param  matrix: matrix to free
  @param  rows: number of rows
  */
  int matrix_i_free (int ***matrix, long rows);

/**
  Writes a double matrix (without parenthesis).
  @param file:       output file
  @param matrix:     matrix to write
  @param rows:       number of rows
  @param columns:    number of columns
  @param tab:        format: leading tabs
  @param separator:  format: separator for columns
  @param ending:     format: end of a row  
  */
  void matrix_d_print (FILE * file, double **matrix, int rows, int columns,
                       char *tab, char *separator, char *ending);

/**
  Writes a double matrix (without parenthesis) with specifically many decimal places.
  @param file:       output file
  @param matrix:     matrix to write
  @param rows:       number of rows
  @param columns:    number of columns
  @param width:      format: number of places altogether
  @param prec:       format: number of decimal places
  @param tab:        format: leading tabs
  @param separator:  format: separator for columns
  @param ending:     format: end of a row
  */
  void matrix_d_print_prec (FILE * file, double **matrix, int rows,
                            int columns, int width, int prec, char *tab,
                            char *separator, char *ending);

/**
  Writes an integer matrix (without parenthesis).
  @param file:       output file
  @param matrix:     matrix to write
  @param rows:       number of rows
  @param columns:    number of columns
  @param tab:        format: leading tabs
  @param separator:  format: separator for columns
  @param ending:     format: end of a row  
  */
  void matrix_i_print (FILE * file, int **matrix, int rows, int columns,
                       char *tab, char *separator, char *ending);

/**
  Reads in a double matrix.
  @return 0 for succes; -1 for error
  @param s:          scanner
  @param matrix:     matrix to read
  @param max_row:    number of rows
  @param max_column: number of columns
  */
  int matrix_d_read (scanner_t * s, double **matrix, int max_row,
                     int max_column);

/**
  Reads in an integer matrix.
  @return 0 for succes; -1 for error
  @param s:          scanner
  @param matrix:     matrix to read
  @param max_row:    number of rows
  @param max_column: number of columns
  */
  int matrix_i_read (scanner_t * s, int **matrix, int max_row,
                     int max_column);

/**
  Determines the number of entries != 0 in a row of a matrix.
  @return         number of entries
  @param matrix:   double matrix
  @param row:      row to scan
  @param max_col:  number of columns 
  */
  int matrix_d_notzero_columns (double **matrix, int row, int max_col);

/**
  Determines the number of entries != 0 in a column of a matrix. 
  @return         number of entries
  @param matrix:   double matrix
  @param col:      column to scan
  @param max_row:  number of rows
  */
  int matrix_d_notzero_rows (double **matrix, int col, int max_row);

/** 
  Scales the rowvectors of a matrix, so that they have sum 1.
  @return 0 for success; -1 for error
  @param matrix:   double matrix
  @param rows:     number of rows
  @param cols:     number of columns
  */
  int matrix_d_normalize (double **matrix, int rows, int cols);

/** 
  Gives the elements in a matrix uniformly distributed values between min and max. 
  @param matrix:   double matrix
  @param rows:     number of rows
  @param cols:     number of columns
  @param min:      minimum for the random values
  @param max:      maximum for the random values
  */
  void matrix_d_random_values (double **matrix, int rows, int cols,
                               double min, double max);

/** 
  Gives the elements in a matrix uniformly distributed values between min and max. 
  Gives all elements in the last row a constant value.
  @param matrix:   double matrix
  @param rows:     number of rows
  @param cols:     number of columns
  @param min:      minimum for the random values
  @param max:      maximum for the random values
  @param c:        value for the last row
  */
  void matrix_d_random_const_values (double **matrix, int rows, int cols,
                                     double min, double max, double c);

/** 
  Gives all elements in a matrix a constant value.
  @param matrix:   double matrix
  @param rows:     number of rows
  @param cols:     number of columns
  @param c:        value for the elements
  */
  void matrix_d_const_values (double **matrix, int rows, int cols, double c);

/** 
  Gives all elements on the 1. upper secondary diagonal in a matrix the value 1.
  @param matrix:   double matrix
  @param rows:     number of rows
  @param cols:     number of columns
  */
  void matrix_d_left_right_strict (double **matrix, int rows, int cols);

/** 
  Gives the elements in a matrix with band width 3 random values. 
  @param matrix:   double matrix
  @param rows:     number of rows
  @param cols:     number of columns
  */
  void matrix_d_random_left_right (double **matrix, int rows, int cols);


/** 
  Gives all elements != 0 in a matrix a constant value.
  @param matrix:   double matrix
  @param rows:     number of rows
  @param cols:     number of columns
  @param c:        value for the elements
  */
  void matrix_d_const_preserve_struct (double **matrix, int rows, int cols,
                                       double c);

/** 
  Gives all elements != 0 in a matrix uniformly distributed random values 
  between 0 and 1.
  @param matrix   double matrix
  @param rows     number of rows
  @param cols     number of columns
  */
  void matrix_d_random_preserve_struct (double **matrix, int rows, int cols);

/** 
  Gives each row in a matrix values according to a certain Gauss density.
  Mean values are randomly generated if mue == NULL.
  u ($\sigma^2$) must be given.
  @return 0 for success; -1 for failure
  @param matrix:   double matrix
  @param rows:     number of rows
  @param cols:     number of columns
  @param mue:      pointer to the vector containing the mean values for each row
  @param u:        standard deviation, for all rows equal
 */
  int matrix_d_gaussrows_values (double **matrix, int rows, int cols,
                                 double **mue, double u);

/** 
  Transposes a matrix.
  @param A:        double matrix
  @param rows:     number of rows
  @param cols:     number of columns
  @param A_T:      transposed double matrix (the returned value)
 */
  void matrix_d_transpose (double **A, int rows, int cols, double **A_T);

/**
  Solves a linear equation system, Ax = b, for a symmetric, positiv definit matrix.
  @return 0 for success; -1 for failure
  @param a:    double matrix 
  @param b:    double vector
  @param dim:  dimension of a
  @param x:    double vector, a solution of the system.
  */
  int matrix_cholesky (double **a, double *b, int dim, double *x);

/**
  Finds the determinant of a symmetric, positiv definit matrix.
  @return 0 for success; -1 for failure
  @param a:    double matrix
  @param dim:  dimension of a
  @param det:  determinant of a, the returning value
  */
  int matrix_det_symposdef (double **a, int dim, double *det);

/** 
  Copies a matrix. Allocation needs to be done outside ! 
  @param src:    double matrix to copy
  @param target: copy of src, the returning value 
  @param rows:   number of rows
  @param cols:   number of columns
*/
  void matrix_d_copy (double **src, double **target, int rows, int cols);

/**
  Checks whether a quadratic double matrix is stochastic
  @return 0/1 flag for true/false
  @param  double NxN matrix to be checked
  @param  matrix dimension N (matrix must be quadaratic)
  */
  int matrix_d_check_stochasticity (double **matrix, int N);



#ifdef __cplusplus
}
#endif
#endif
/*@} (Doc++-Group: matrix) */
