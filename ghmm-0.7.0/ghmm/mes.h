/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/mes.h
*       Authors:  Frank Nuebel, Benjamin Georgi
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
*       This file is version $Revision: 1333 $ 
*                       from $Date: 2005-09-12 23:39:46 +0200 (Mon, 12 Sep 2005) $
*             last change by $Author: grunau $.
*
*******************************************************************************/
#ifndef MES_H
#define MES_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ghmm/internal.h>

#ifdef __cplusplus
extern "C" {
#endif /*__cplusplus*/

/**
   @name Help functions concerning memory and small calculations.
*/

#define MES_0_PTR 0
#define MES_NEG_ARG 1
#define MES_BIG_ARG 2
#define MES_0_ARG 3
#define MES_FALSE_ARG 4

#define MES_FLAG_TIME        1
#define MES_FLAG_WIN         4
#define MES_FLAG_FILE       16
#define MES_FLAG_FILE_WIN   ( MES_FLAG_FILE | MES_FLAG_WIN )
#define MES_FLAG_TIME_WIN   ( MES_FLAG_TIME | MES_FLAG_FILE | MES_FLAG_WIN )

#define MES_MODUL_INFO  "(" __DATE__ ":" __FILE__ ")"
#define MES_PROC_INFO   "(" __DATE__ ":" __FILE__  ":" CUR_PROC ")"


#define MES_WIN             MES_FLAG_WIN, -1, NULL, NULL
#define MES_FILE            MES_FLAG_FILE, -1, NULL, NULL
#define MES_TIME            (MES_FLAG_TIME | MES_FLAG_FILE), -1, NULL, NULL
#define MES_FILE_WIN        MES_FLAG_FILE_WIN, -1, NULL, NULL
#define MES_TIME_WIN        MES_FLAG_TIME_WIN, -1, NULL, NULL
#define MES_PROT_WIN        MES_FLAG_WIN, __LINE__ ,MES_PROC_INFO, CUR_PROC
#define MES_PROT_FILE       MES_FLAG_FILE, __LINE__ ,MES_PROC_INFO, CUR_PROC
#define MES_PROT            MES_FLAG_FILE_WIN, __LINE__ ,MES_PROC_INFO, CUR_PROC
#define MES_PROT_TIME       MES_FLAG_TIME_WIN, __LINE__ ,MES_PROC_INFO, CUR_PROC

/* stuff from sys.h */

#if defined(_WIN32)
#  define getthreadid() GetCurrentThreadId()
#  define getprocessid() GetCurrentProcessId()
#else
#  include <unistd.h>
#  define getthreadid()  1
#  define getprocessid() 1
#endif                          /* defined(WIN32) */


/* end stuff from sys.h */

/* stuff from stdmacro.h */

#ifndef strtoupper
#define strtoupper(str) if(1){char*p=(str);if(p)for(;*p;p++)*p=toupper(*p);}else
#endif

#ifndef strtolower
#define strtolower(str) if(1){char*p=(str);if(p)for(;*p;p++)*p=tolower(*p);}else
#endif

#ifndef m_max
#define m_max( a, b ) ( ( (a) > (b) ) ? (a) : (b) )
#endif

#ifndef m_min
#define m_min( a, b ) ( ( (a) < (b) ) ? (a) : (b) )
#endif

#ifndef m_int
#define m_int( x ) ( ((x)>=0) ? (int)((x)+0.5) : (int)((x)-0.5) )
#endif

#ifndef m_integer
#define m_integer( x ) ( ((x)>=0) ? floor(x) : ceil(x) )
#endif

#ifndef m_frac
#define m_frac( x ) ( (x) - m_integer(x) )
#endif

#ifndef m_sqr
#define m_sqr( x ) ( (x)*(x) )
#endif

#ifndef m_bitget
#define m_bitget( BitStr, Bit ) (((char*)(BitStr))[(Bit)>>3] &  (1<<((Bit)&7)) )
#endif

#ifndef m_bitset
#define m_bitset( BitStr, Bit ) (((char*)(BitStr))[(Bit)>>3] |= (1<<((Bit)&7)) )
#endif

#ifndef m_bitclr
#define m_bitclr( BitStr, Bit ) (((char*)(BitStr))[(Bit)>>3] &= ~(1<<((Bit)&7)))
#endif

#ifndef m_bitinv
#define m_bitinv( BitStr, Bit ) (((char*)(BitStr))[(Bit)>>3] ^= ~(1<<((Bit)&7)))
#endif

/* */
#ifndef m_free
#define m_free( p )  {if(p) { free(p); (p) = NULL; } else { printf("ERROR: Attempted m_free on NULL pointer.  Bad program. BAD ! No cookie for you.\n\n");abort();}}
#endif

#ifndef m_strlen
#define m_strlen( src, cnt ) ( (src) ? ( (cnt<0) ? strlen(src) : cnt) : -1 )
#endif

#ifndef m_align
#define m_align( a, n ) ( ((n)&&((a)%(n))) ? ((a)+(n)-((a)%(n))) : (a) )
#endif

#ifndef m_ispow2
#define m_ispow2( n ) ( !((n)&((n)-1)) )
#endif

#ifndef m_gauss
#define m_gauss( a, n ) ( ((n)&&((a)%(n))) ? ((a)-((a)%(n))) : (a) )
#endif

#ifndef m_approx
#define m_approx( a, b, eps ) ( ( (a) - (eps) <=  (b) && (a) + (eps) >= (b) ) ? 1 : 0)
#endif

#ifndef m_ptr
#define m_ptr( p, offs ) ( (void*)((char*)(p)+(offs)) )
#endif

#ifndef m_memset
#define m_memset(dst, c, entries) memset( (dst), (c), sizeof(*(dst))*(entries))
#endif

#ifndef m_memcpy
#define m_memcpy(dst, src, entries) memcpy((dst),(src),sizeof(*(dst))*(entries))
#endif

#ifndef m_fclose
#define m_fclose( fp ) \
  if((fp) && (fp) - stdout) { fclose(fp); (fp) = NULL; } else
#endif

#ifndef m_fread
#define m_fread( fp, dest, cnt )  mes_fread((fp),(dest),(cnt)*sizeof(*(dest)) )
#endif

#ifndef m_fwrite
#define m_fwrite( fp, src, cnt )  mes_fwrite((fp),(src),(cnt)*sizeof(*(src)) )
#endif

/* neu fuer hmm: ungefaehr gleiche Werte, (BW) */
#ifndef m_approx
#define m_approx( a, b, eps ) ( ( (a) - (eps) <=  (b) && (a) + (eps) >= (b) ) ? 1 : 0)
#endif

/* end of things from stdmacro.h */

  /**
   */
#define mes_check_ptr(arg, call) \
  if (!(arg)) { mes_err( #arg, MES_0_PTR, MES_PROC_INFO ); call; }
  /**
   */
#define mes_check_bnd( arg, bnd, call ) \
  if((arg) > (bnd)){ mes_err( #arg, MES_BIG_ARG, MES_PROC_INFO ); call;} else
  /**
   */
#define mes_check_0( arg, call ) \
  if(!(arg)) { mes_err( #arg, MES_0_ARG, MES_PROC_INFO ); call;} else
  /**
   */
#define mes_check_neg( arg, call ) \
  if((arg)<0){ mes_err( #arg, MES_NEG_ARG, MES_PROC_INFO );call;} else
  /**
   */
#define mes_check_expr( arg, call ) \
  if(!(arg)){ mes_err( #arg, MES_FALSE_ARG, MES_PROC_INFO );call;} else


  /**
   */
#define mes_file( txt )     mes_smart( MES_FLAG_FILE, txt, -1 )
  /**
   */
#define mes_file_win( txt ) mes_smart( MES_FLAG_FILE_WIN, txt, -1 )
  /**
   */
#define mes_win( txt )      mes_smart( MES_FLAG_WIN, txt, -1 )

  /**
   */
#define mes_proc()      mes( MES_PROT, NULL )
  /**
   */
#define mes_prot( txt ) mes( MES_PROT_TIME, txt )

  /**
   */
  char *mes_get_std_path (void);

  /**
   */
  void mes (int flags, int line, char *xproc, char *proc, char *format, ...);
  /**
   */
  void mes_smart (int flags, const char *txt, int bytes);
  /**
   */
  int mes_ability (int on);
  /**
   */
  void mes_err (char *txt, int error_nr, char *proc_info);
  /**
   */
  void mes_exit (void);
  /**
   */
  void mes_fformat (char *txt, char *filename, int line, char *proc_info);
  /**
   */
  void mes_init (char *logfile, void (*winfct) (const char *), int argc,
                 char *argv[]);
  /**
   */
  void mes_init_args (int argc, char *argv[]);
  /**
   */
  void mes_init_logfile (char *file_name);
  /**
   */
  void mes_init_winfct (void (*win_fct) (const char *));
  /**
   */
  int mes_insert (FILE * fp, char src, int cnt);
  /**
   */
  void mes_time (void);
  /**
   */
  void mes_va (int flags, int line, char *xproc, char *proc, char *format,
               va_list args);
  /**
   */
  int mes_win_ability (int on);

  /**
   */
  void *mes_calloc (int bytes);
  /**
   */
  int mes_copy (char *oldname, char *newname);
  /**
   */
  int mes_fgetc (FILE * fp);
  /**
   */
  int mes_fflush (FILE * fp);
  /**
   */
  FILE *mes_fopen (const char *filename, char *attribute_string);
  /**
   */
  int mes_fprintf (FILE * fp, char *format, ...);
  /**
   */
  int mes_fputc (FILE * fp, char chr);
  /**
   */
  int mes_fputs (FILE * fp, char *str);
  /**
   */
  int mes_fread (FILE * fp, void *mem, int bytes);
  /**
   */
  int mes_fread_quiet (FILE * fp, void *mem, int bytes);
  /**
   */
  int mes_fseek (FILE * fp, long offset, int fromwhere);
#ifdef WIN32
  /**
   */
  int mes_fseek64 (FILE * fp, unsigned int uoff, unsigned int loff,
                   int fromwhere);
#endif
  /**
   */
  int mes_ftell (FILE * fp);
  /**
   */
  int mes_fwrite (FILE * fp, void *mem, int bytes);
  /**
   */
  void *mes_malloc (int bytes);
  /**
   */
  int mes_move (char *oldname, char *newname);
  /**
   */
  int mes_realloc (void **mem, int bytes);
  /**
   */
  int mes_remove (char *filename);
  /**
   */
  int mes_rename (char *oldname, char *newname);
  /**
   */
  FILE *mes_tmpfile (void);

#ifdef __cplusplus
}
#endif /*__cplusplus*/
#endif                          /* MES_H */
