/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/gauss_tail.c
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
*       This file is version $Revision: 1198 $
*                       from $Date: 2005-07-05 17:38:09 +0200 (Tue, 05 Jul 2005) $
*             last change by $Author: schliep $.
*
*******************************************************************************/
   
  
#include <float.h>
#include <math.h>
#include "const.h"
#include "randvar.h"
  
/*============================================================================*/ 
double pmue (double mue, double A, double B, double eps)
{
  
double feps, u, Atil, Btil;
  
Atil = A + eps;
  
Btil = B + eps * A;
  
u = Btil - mue * Atil;
  
    /* if (u < EPS_U) u = (double)EPS_U; DANGEROUS: would fudge the function value! */ 
    if (u <= DBL_MIN)
    
return (mue - A);
  
feps = randvar_normal_density_trunc (-eps, mue, u, -eps);
  
return (A - mue - u * feps);

}



/*============================================================================*/ 
/* To avoid numerical ocillation:
   Interpolate p(\mu) between 2 sampling points for PHI
   NOTA BENE: This Version is very expensive and exact. 
*/ 
double pmue_interpol (double mue, double A, double B, double eps)
{
  
double u, Atil, Btil, z, z1, z2, m1, m2, u1, u2, p1, p2, pz;
  
int i1, i2;
  
Atil = A + eps;
  
Btil = B + eps * A;
  
u = Btil - mue * Atil;
  
    /*if (u < EPS_U) u = (double)EPS_U; DANGEROUS: would fudge the function value! */ 
    if (u <= DBL_MIN)
    
return (mue - A);
  

    /* Compute like normally where mue positiv. */ 
    if (mue >= 0.0)
    
return (A - mue - u * randvar_normal_density_trunc (-eps, mue, u, -eps));
  

    /* Otherwise: Interpolate the function itself between 2 sampling points. */ 
    z = (eps + mue) / sqrt (u);
  

i1 = (int) (fabs (z) * randvar_get_xfaktphi ());
  
if (i1 >= randvar_get_philen () - 1) {
    
i1 = (int) randvar_get_philen () - 1;
    
i2 = i1;
  
}
  
  else
    
i2 = i1 + 1;
  
z1 = i1 / randvar_get_xfaktphi ();
  
z2 = i2 / randvar_get_xfaktphi ();
  

m1 = -z1 * sqrt (Btil + eps * Atil + Atil * Atil * z1 * z1 * 0.25) 
    -(eps + Atil * z1 * z1 * 0.5);
  
m2 = -z2 * sqrt (Btil + eps * Atil + Atil * Atil * z2 * z2 * 0.25) 
    -(eps + Atil * z2 * z2 * 0.5);
  
u1 = Btil - m1 * Atil;
  
u2 = Btil - m2 * Atil;
  

p1 = A - m1 - u1 * randvar_normal_density_trunc (-eps, m1, u1, -eps);
  
p2 = A - m1 - u1 * randvar_normal_density_trunc (-eps, m2, u2, -eps);
  

if (i1 >= randvar_get_philen () - 1)
    
pz = p1;
  
  else {
    
pz = p1 + (fabs (z) - i1 * randvar_get_xstepphi ()) 
      *(p2 - p1) / randvar_get_xstepphi ();
    
      /* pz = p1; */ 
  }
  
return (pz);

}



/*============================================================================*/ 
double pmue_umin (double mue, double A, double B, double eps)
{
  
double feps, u;
  
u = EPS_U;
  
feps = randvar_normal_density_trunc (-eps, mue, u, -eps);
  
return (A - mue - u * feps);

}


