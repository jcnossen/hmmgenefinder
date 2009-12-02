/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/randvar.c
*       Authors:  Bernhard Knab, Benjamin Rich
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
#include <string.h>
#ifdef HAVE_LIBPTHREAD
# include <pthread.h>
#endif /* HAVE_LIBPTHREAD */
#include "mes.h"
#include "mprintf.h"
#include "randvar.h"
#include "rng.h"
#include "float.h"
#include "const.h"

#ifdef DO_WITH_GSL

# include <gsl/gsl_math.h>
# include <gsl/gsl_sf_erf.h>
# include <gsl/gsl_randist.h>

/* missing functions in gsl-0.7 */

#define EVAL_RESULT(fn) \
   gsl_sf_result result; \
   int status = fn; \
   if (status == GSL_EDOM) { \
     return GSL_NAN; \
   } else if (status != GSL_SUCCESS) { \
     GSL_ERROR(#fn, status); \
   } ; \
   return result.val;

# ifndef HAVE_GSL_SF_ERF
double gsl_sf_erfc (double x)
{
  EVAL_RESULT (gsl_sf_erfc_e (x, &result));
}
# endif

# ifndef HAVE_GSL_SF_ERFC
double gsl_sf_erf (double x)
{
  EVAL_RESULT (gsl_sf_erf_e (x, &result));
}
# endif

#else

#include "scanner.h"

#endif /* DO_WITH_GSL */

/* A list of already calculated values of the density function of a 
   N(0,1)-distribution, with x in [0.00, 19.99] */
#define PDFLEN 2000
#define X_STEP_PDF 0.01         /* step size */
#define X_FAKT_PDF 100          /* equivalent to step size */
static double pdf_stdnormal[PDFLEN];
static int pdf_stdnormal_exists = 0;

/* A list of already calulated values PHI of the Gauss distribution is
   read in, x in [-9.999, 0] */
#define X_STEP_PHI 0.001        /* step size */
#define X_FAKT_PHI 1000         /* equivalent to step size */
static int PHI_len = 0;
#if 0                           /* needed for randvar_get_xPHIxgleichPHIy() */
static double x_PHI_xy = -1.0;
#endif
static double x_PHI_1 = -1.0;


#ifndef DO_WITH_GSL
static double *PHI = NULL;


static int randvar_read_PHI ()
{
# define CUR_PROC "randvar_read_PHI"
  int res = -1;
  char filename = "ghmm/PHI_001_20.dat";
  scanner_t *s = NULL;

//#warning "PHI_DATA_FILE deprecated!"
  s = scanner_alloc (filename);
  if (!s) {
    mes_proc ();
    goto STOP;
  }
  scanner_get_name (s);
  scanner_consume (s, '=');
  if (s->err)
    goto STOP;
  if (!strcmp (s->id, "PHI")) {
    scanner_consume (s, '{');
    if (s->err)
      goto STOP;
    PHI = scanner_get_double_earray (s, &PHI_len);
    if (s->err)
      goto STOP;
    scanner_consume (s, ';');
    if (s->err)
      goto STOP;
    scanner_consume (s, '}');
    if (s->err)
      goto STOP;
    scanner_consume (s, ';');
    if (s->err)
      goto STOP;
  }
  else {
    scanner_error (s, "unknown identifier");
    goto STOP;
  }
  // printf("%.4f\n", PHI[PHI_len-1]); 

  res = 0;;
STOP:
  scanner_free (&s);
  return res;
# undef CUR_PROC
}       



/* randvar_read_PHI */


/*----------------------------------------------------------------------------*/
static int randvar_init_PHI ()
{
# define CUR_PROC "randvar_init_PHI"
#ifdef HAVE_LIBPTHREAD
  static pthread_mutex_t lock;
#endif /* HAVE_LIBPTHREAD */
  /* Read PHI in */
  if (!PHI_len) {
#ifdef HAVE_LIBPTHREAD
    pthread_mutex_lock (&lock); /* Put on a lock, because the clustering is parallel */
#endif /* HAVE_LIBPTHREAD */
    if (randvar_read_PHI () == -1) {
      mes_proc ();
      goto STOP;
    };
#ifdef HAVE_LIBPTHREAD
    pthread_mutex_unlock (&lock);       /* Take the lock off */
#endif /* HAVE_LIBPTHREAD */
  }
  return 0;
STOP:
  return (-1);
# undef CUR_PROC
}                               /* randvar_init_PHI */


#endif /* DO_WITH_GSL */


/*============================================================================*/
/* needed by pmue_interpol */

double randvar_get_xfaktphi ()
{
  return X_FAKT_PHI;
}

double randvar_get_xstepphi ()
{
  return X_STEP_PHI;
}

double randvar_get_philen ()
{
#ifdef DO_WITH_GSL
  return PHI_len;
#else
  return randvar_get_xPHIless1 () / X_STEP_PHI;
#endif
}


/*============================================================================*/
double randvar_get_PHI (double x)
{
# define CUR_PROC "randvar_get_PHI"

#ifdef DO_WITH_GSL
  return (gsl_sf_erf (x * M_SQRT1_2) + 1.0) / 2.0;
#else
  int i;
  double phi_x;

  if (randvar_init_PHI () == -1) {
    mes_proc ();
    goto STOP;
  }

  /* Linear interpolation (Alternative: Round off with i=m_int(fabs(x)*X_FAKT)) */
  i = (int) (fabs (x) * X_FAKT_PHI);
  if (i >= PHI_len - 1) {
    i = PHI_len - 1;
    phi_x = PHI[i];
  }
  else
    phi_x =
      PHI[i] + (fabs (x) - i * X_STEP_PHI) * (PHI[i + 1] -
                                              PHI[i]) / X_STEP_PHI;
  /* NOTA BENE: PHI is tabulated for negative values! */
  if (x > 0.0)
    return (1.0 - phi_x);
  else
    return (phi_x);
#endif /* DO_WITH_GSL */
STOP:
  return (-1.0);
# undef CUR_PROC
}                               /* randvar_get_PHI */


/*============================================================================*/
/* When is PHI[x,0,1] == 1? */
double randvar_get_xPHIless1 ()
{
# define CUR_PROC "randvar_get_xPHIless1"
#ifdef DO_WITH_GSL
  if (x_PHI_1 == -1) {
    double low, up, half;
    low = 0;
    up = 100;
    while (up - low > 0.001) {
      half = (low + up) / 2.0;
      if (randvar_get_PHI (half) < 1.0)
        low = half;
      else
        up = half;
    }
    x_PHI_1 = low;
  }
  return (x_PHI_1);
#else
  double x;
  int i;
  if (x_PHI_1 == -1) {
    if (randvar_init_PHI () == -1) {
      mes_proc ();
      goto STOP;
    }
    /* The last value of the table is 1 */
    for (x = (PHI_len - 1) * X_STEP_PHI, i = PHI_len - 1; i > 0;
         x -= X_STEP_PHI, i--)
      if (randvar_get_PHI (-x) > 0.0)
        break;
    /* Modification: x exactly between 2 sampling points! */
    x_PHI_1 = x - (double) X_STEP_PHI / 2.0;
  }
  return (x_PHI_1);
#endif
STOP:
  return (-1.0);
# undef CUR_PROC
}

#if 0                           /* is not in use */
/*============================================================================*/
/* When is PHI[x,0,1] ==  PHI[y,0,1] for consecutive x and y ?*/
double randvar_get_xPHIxgleichPHIy ()
{
# define CUR_PROC "randvar_get_xPHIxgleichPHIy"
  double x, y;
  int i;
  if (x_PHI_xy == -1) {
    if (randvar_init_PHI () == -1) {
      mes_proc ();
      goto STOP;
    }
    y = -1.0;
    for (x = 0.0, i = 0; i < PHI_len; x += X_STEP_PHI, i++) {
      if (randvar_get_PHI (-x) == randvar_get_PHI (-y))
        break;
      y = x;
    }
    x_PHI_xy = y;
  }
  return (x_PHI_xy);
STOP:
  return (-1.0);
# undef CUR_PROC
}
#endif /* 0 */

/*============================================================================*/
double randvar_get_1overa (double x, double mean, double u)
{
  /* Calulates 1/a(x, mean, u), with a = the integral from x til \infty over
     the Gauss density function */
# define CUR_PROC "randvar_get_1overa"

#ifdef DO_WITH_GSL
  double erfc_value;
#else
  int i;
  double y, z, phi_z, a;
#endif

  if (u <= 0.0) {
    mes_prot ("u <= 0.0 not allowed\n");
    goto STOP;
  }

#ifdef DO_WITH_GSL
  /* int gsl_sf_erfc (double x) 
     erfc(x) = 1 - erf(x) = 2/\sqrt(\pi) \int_x^\infty \exp(-t^2)
   */
  erfc_value = gsl_sf_erfc ((x - mean) / sqrt (u * 2));
  if (erfc_value <= DBL_MIN) {
    mes (MES_WIN, "a ~= 0.0 critical! (mue = %.2f, u =%.2f)\n", mean, u);
    return (erfc_value);
  }
  else
    return (2.0 / erfc_value);
#else

  if (randvar_init_PHI () == -1) {
    mes_proc ();
    goto STOP;
  };

  y = 1 / sqrt (u);
  z = (x - mean) * y;
  /* Linear interpolation (Alternative: Round off with i=m_int(fabs(z)*X_FAKT)) */
  i = (int) (fabs (z) * X_FAKT_PHI);

  if (i >= PHI_len - 1) {
    i = PHI_len - 2;
    /* Originally:
       i = PHI_len-1; but then, the last value in the table is zero! */
    phi_z = PHI[i];
  }
  else
    phi_z =
      PHI[i] + (fabs (z) - i * X_STEP_PHI) * (PHI[i + 1] -
                                              PHI[i]) / X_STEP_PHI;
  /* NOTA BENE: PHI is tabulated for negative values! */
  if (z > 0.0) {
    if (phi_z == 0) {
      mes_proc ();
      goto STOP;
    }
    else
      a = 1 / phi_z;            /* PHI between 0.5 and 1 */ /* ??? between 0.5 and 0 ! */
  }
  else {
    a = 1 - phi_z;
    if (a > DBL_MIN)
      a = 1 / a;
    else {
      a = 0.0;
      mes (MES_WIN, "a ~= 0.0 critical! (mue = %.2f, u =%.2f)\n", mean, u);     /* goto STOP; */
    }
  }
  return a;
#endif

STOP:
  return (-1.0);
# undef CUR_PROC
}                               /* randvar_get_1overa */


/*============================================================================*/
/* REMARK:
   The calulation of this density function was testet, by calculating the 
   following integral sum for arbitrary mue and u:
     for (x = 0, x < ..., x += step(=0.01/0.001/0.0001)) 
       isum += step * randvar_normal_density_pos(x, mue, u);
   In each case, the sum "converged" evidently towards 1!
   (BK, 14.6.99)
   CHANGE:
   Truncate at -EPS_NDT (const.h), so that x = 0 doesn't lead to a problem.
   (BK, 15.3.2000)
*/
double randvar_normal_density_pos (double x, double mean, double u)
{
# define CUR_PROC "randvar_normal_density_pos"
  return randvar_normal_density_trunc (x, mean, u, -EPS_NDT);
# undef CUR_PROC
}                               /* double randvar_normal_density_pos */


/*============================================================================*/
double randvar_normal_density_trunc (double x, double mean, double u,
                                     double a)
{
# define CUR_PROC "randvar_normal_density_trunc"
#ifndef DO_WITH_GSL
  double c;
#endif /* DO_WITH_GSL */

  if (u <= 0.0) {
    mes_prot ("u <= 0.0 not allowed\n");
    goto STOP;
  }
  if (x < a)
    return (0.0);

#ifdef DO_WITH_GSL
  /* move mean to the right position */
  /* double gsl_ran_gaussian_tail_pdf (double x, double a, double sigma) */
  return gsl_ran_gaussian_tail_pdf (x - mean, a - mean, sqrt (u));
#else
  if ((c = randvar_get_1overa (a, mean, u)) == -1) {
    mes_proc ();
    goto STOP;
  };
  return (c * randvar_normal_density (x, mean, u));
#endif /* DO_WITH_GSL */

STOP:
  return (-1.0);
# undef CUR_PROC
}                               /* double randvar_normal_density_trunc */


/*============================================================================*/
double randvar_normal_density (double x, double mean, double u)
{
# define CUR_PROC "randvar_normal_density"
#ifndef DO_WITH_GSL
  double expo;
#endif
  if (u <= 0.0) {
    mes_prot ("u <= 0.0 not allowed\n");
    goto STOP;
  }
  /* The denominator is possibly < EPS??? Check that ? */
#ifdef DO_WITH_GSL
  /* double gsl_ran_gaussian_pdf (double x, double sigma) */
  return gsl_ran_gaussian_pdf (x - mean, sqrt (u));
#else
  expo = exp (-1 * m_sqr (mean - x) / (2 * u));
  return (1 / (sqrt (2 * PI * u)) * expo);
#endif

STOP:
  return (-1.0);
# undef CUR_PROC
}                               /* double randvar_normal_density */


/*============================================================================*/
/* special smodel pdf need it: smo->density==normal_approx: */
/* generates a table of of aequidistant samples of gaussian pdf */

static int randvar_init_pdf_stdnormal ()
{
# define CUR_PROC "randvar_init_pdf_stdnormal"
  int i;
  double x = 0.00;
  for (i = 0; i < PDFLEN; i++) {
    pdf_stdnormal[i] = 1 / (sqrt (2 * PI)) * exp (-1 * x * x / 2);
    x += (double) X_STEP_PDF;
  }
  pdf_stdnormal_exists = 1;
  /* printf("pdf_stdnormal_exists = %d\n", pdf_stdnormal_exists); */
  return (0);
# undef CUR_PROC
}                               /* randvar_init_pdf_stdnormal */


double randvar_normal_density_approx (double x, double mean, double u)
{
# define CUR_PROC "randvar_normal_density_approx"
#ifdef HAVE_LIBPTHREAD
  static pthread_mutex_t lock;
#endif /* HAVE_LIBPTHREAD */
  int i;
  double y, z, pdf_x;
  if (u <= 0.0) {
    mes_prot ("u <= 0.0 not allowed\n");
    goto STOP;
  }
  if (!pdf_stdnormal_exists) {
#ifdef HAVE_LIBPTHREAD
    pthread_mutex_lock (&lock); /* Put on a lock, because the clustering is parallel   */
#endif /* HAVE_LIBPTHREAD */
    randvar_init_pdf_stdnormal ();
#ifdef HAVE_LIBPTHREAD
    pthread_mutex_unlock (&lock);       /* Take the lock off */
#endif /* HAVE_LIBPTHREAD */
  }
  y = 1 / sqrt (u);
  z = fabs ((x - mean) * y);
  i = (int) (z * X_FAKT_PDF);
  /* linear interpolation: */
  if (i >= PDFLEN - 1) {
    i = PDFLEN - 1;
    pdf_x = y * pdf_stdnormal[i];
  }
  else
    pdf_x = y * (pdf_stdnormal[i] +
                 (z - i * X_STEP_PDF) *
                 (pdf_stdnormal[i + 1] - pdf_stdnormal[i]) / X_STEP_PDF);
  return (pdf_x);
STOP:
  return (-1.0);
# undef CUR_PROC
}                               /* double randvar_normal_density_approx */


/*============================================================================*/
double randvar_std_normal (int seed)
{
# define CUR_PROC "randvar_std_normal"
  if (seed != 0) {
    GHMM_RNG_SET (RNG, seed);
    return (1.0);
  }
  else {
#ifdef DO_WITH_GSL
    return (gsl_ran_gaussian (RNG, 1.0));
#else
    /* Use the polar Box-Mueller transform */
    /*
       double x, y, r2;

       do {
       x = 2.0 * GHMM_RNG_UNIFORM(RNG) - 1.0;
       y = 2.0 * GHMM_RNG_UNIFORM(RNG) - 1.0;
       r2 = (x * x) + (y * y);
       } while (r2 >= 1.0);

       return x * sqrt((-2.0 * log(r2)) / r2);
     */

    double r2, theta;

    r2 = -2.0 * log (GHMM_RNG_UNIFORM (RNG));   /* r2 ~ chi-square(2) */
    theta = 2.0 * PI * GHMM_RNG_UNIFORM (RNG);  /* theta ~ uniform(0, 2 \pi) */
    return sqrt (r2) * cos (theta);
#endif
  }
# undef CUR_PROC
}                               /* randvar_std_normal */


/*============================================================================*/
double randvar_normal (double mue, double u, int seed)
{
# define CUR_PROC "randvar_normal"
  if (seed != 0) {
    GHMM_RNG_SET (RNG, seed);
    return (1.0 * sqrt (u) + mue);
  }
  else {
#ifdef DO_WITH_GSL
    return (gsl_ran_gaussian (RNG, sqrt (u)) + mue);
#else
    double x;
    x = sqrt (u) * randvar_std_normal (seed) + mue;
    return (x);
#endif
  }
# undef CUR_PROC
}                               /* randvar_normal */


/*============================================================================*/
#ifndef DO_WITH_GSL
# define C0 2.515517
# define C1 0.802853
# define C2 0.010328
# define D1 1.432788
# define D2 0.189269
# define D3 0.001308
#endif

double randvar_normal_pos (double mue, double u, int seed)
{
# define CUR_PROC "randvar_normal_pos"
  double x = -1;
  double sigma;
#ifdef DO_WITH_GSL
  double s;
#else
  double U, Us, Us1, Feps, Feps1, t, T;
#endif

  if (u <= 0.0) {
    mes_prot ("u <= 0.0 not allowed\n");
    goto STOP;
  }
  sigma = sqrt (u);

  if (seed != 0) {
    GHMM_RNG_SET (RNG, seed);
    return (1.0);
  }

#ifdef DO_WITH_GSL
  /* up to version 0.8 gsl_ran_gaussian_tail can not handle negative cutoff */
#define GSL_RAN_GAUSSIAN_TAIL_BUG 1
#ifdef GSL_RAN_GAUSSIAN_TAIL_BUG
  s = (-mue) / sigma;
  if (s < 1) {
    do {
      x = gsl_ran_gaussian (RNG, 1.0);
    }
    while (x < s);
    return x * sigma + mue;
  }
#endif /* GSL_RAN_GAUSSIAN_TAIL_BUG */
  /* move boundary to lower values in order to achieve maximum at mue
     gsl_ran_gaussian_tail(generator, lower_boundary, sigma)
   */
  return gsl_ran_gaussian_tail (RNG, -mue, sqrt (u)) + mue;

#else /* DO_WITH_GSL */
  /* Method: Generate Gauss-distributed random nunbers (with GSL-lib.),
     until a positive one is found -> not very effective if mue << 0
     while (x < 0.0) {
     x = sigma * randvar_std_normal(seed) + mue;
     } */

  /* Inverse transformation with restricted sampling by Fishman */
  U = GHMM_RNG_UNIFORM (RNG);
  Feps = randvar_get_PHI (-(EPS_NDT + mue) / sigma);
  Us = Feps + (1 - Feps) * U;
  /* Numerically better: 1-Us = 1-Feps - (1-Feps)*U, therefore: 
     Feps1 = 1-Feps, Us1 = 1-Us */
  Feps1 = randvar_get_PHI ((EPS_NDT + mue) / sigma);
  Us1 = Feps1 - Feps1 * U;
  t = m_min (Us, Us1);
  t = sqrt (-log (t * t));
  T =
    sigma * (t -
             (C0 + t * (C1 + t * C2)) / (1 + t * (D1 + t * (D2 + t * D3))));
  if (Us - 0.5 < 0)
    x = mue - T;
  else
    x = mue + T;
#endif /* DO_WITH_GSL */


STOP:
  return (x);
# undef CUR_PROC
}                               /* randvar_normal_pos */


/*============================================================================*/
double randvar_uniform_int (int seed, int K)
{
# define CUR_PROC "randvar_uniform_int"
  if (seed != 0) {
    GHMM_RNG_SET (RNG, seed);
    return (1.0);
  }
  else {
#ifdef DO_WITH_GSL
    /* more direct solution than old version ! */
    return (double) gsl_rng_uniform_int (RNG, K);
#else
    return (double) ((int) (((double) K) * GHMM_RNG_UNIFORM (RNG)));
#endif
  }
# undef CUR_PROC
}                               /* randvar_uniform_int */


/*============================================================================*/
/* cumalative distribution function of N(mean, u) */
double randvar_normal_cdf (double x, double mean, double u)
{
# define CUR_PROC "randvar_normal_cdf"
  if (u <= 0.0) {
    mes_prot ("u <= 0.0 not allowed\n");
    goto STOP;
  }
#ifdef DO_WITH_GSL
  /* PHI(x)=erf(x/sqrt(2))/2+0.5 */
  return (gsl_sf_erf ((x - mean) / sqrt (u * 2.0)) + 1.0) / 2.0;
#else
  /* The denominator is possibly < EPS??? Check that ? */
  return (randvar_get_PHI ((x - mean) / sqrt (u)));
#endif /* DO_WITH_GSL */
STOP:
  return (-1.0);
# undef CUR_PROC
}                               /* double randvar_normal_cdf */

/*============================================================================*/
/* cumalative distribution function of -EPS_NDT-truncated N(mean, u) */
double randvar_normal_pos_cdf (double x, double mean, double u)
{
# define CUR_PROC "randvar_normal_pos_cdf"
#ifndef DO_WITH_GSL
  double Fx, c;
#endif

  if (x <= 0.0)
    return (0.0);
  if (u <= 0.0) {
    mes_prot ("u <= 0.0 not allowed\n");
    goto STOP;
  }
#ifdef DO_WITH_GSL
  /*
     Function: int gsl_sf_erfc_e (double x, gsl_sf_result * result) 
     These routines compute the complementary error function
     erfc(x) = 1 - erf(x) = 2/\sqrt(\pi) \int_x^\infty \exp(-t^2). 
   */
  return 1.0 + (gsl_sf_erf ((x - mean) / sqrt (u * 2)) -
                1.0) / gsl_sf_erfc ((-mean) / sqrt (u * 2));
#else
  /*The denominator is possibly < EPS??? Check that ? */
  Fx = randvar_get_PHI ((x - mean) / sqrt (u));
  c = randvar_get_1overa (-EPS_NDT, mean, u);
  return (c * (Fx - 1) + 1);
#endif /* DO_WITH_GSL */
STOP:
  return (-1.0);
# undef CUR_PROC
}                               /* double randvar_normal_cdf */
