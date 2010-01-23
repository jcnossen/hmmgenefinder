/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/mes.c
*       Authors:  Frank Nübel, Benjamin Georgi
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
*       This file is version $Revision: 1333 $
*                       from $Date: 2005-09-12 23:39:46 +0200 (Mon, 12 Sep 2005) $
*             last change by $Author: grunau $.
*
*******************************************************************************/

#include <errno.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#  include <windows.h>
#  include <io.h>
#endif

#include "mes.h"
#include "mprintf.h"

/* from sys.h */

#if defined( _WINDOWS ) || defined( _DOS ) || defined(_PPC)
#  define strcasecmp( a, b ) strcmpi( a, b )
#  define PathSepChar ('\\')
#else
#  define PathSepChar ('/')
#endif


#define MAX_MES_FILE_LENGTH 1400000

static char *mes_err_txt[] = {
  "0 pointer",
  "negative value",
  "value out of range",
  "zero value",
  "boolean value FALSE"
};

typedef struct mes_t {
  int thread_id;
  char *std_path;
  char *logfile;
  int argc;
  char **argv;
  char win_enabled;
  char enabled;
  void (*win_fkt) (const char *);
} mes_t;

#define MES_PROCESS_MAX 0x10000
static mes_t *mes_process[MES_PROCESS_MAX];
static int mes_process_n = 0;

static char mes_default_logfile[] = "message.txt";
static char mes_default_bak_logfile[] = "message.bak";

/*static char*  mes_logfile = NULL;
static int    mes_argc = 0;
static char** mes_argv = NULL;

static char  mes_win_enabled = 1;
static char  mes_enabled = 1;
static void(*mes_win_fkt)(char*) = NULL;*/

/*---------------------------------------------------------------------------*/

static mes_t *mes_process_get (void)
{
  int i, tid = getthreadid ();
  for (i = mes_process_n - 1; i >= 0; i--)
    if (mes_process[i])
      if (mes_process[i]->thread_id == tid)
        return (mes_process[i]);
  return (NULL);
}                               /* mes_process_get */

/*---------------------------------------------------------------------------*/
static int mes_process_enabled (void)
{
  mes_t *mpc = mes_process_get ();
  if (mpc)
    return mpc->enabled;
  return (1);
}                               /* mes_process_enabled */

/*---------------------------------------------------------------------------*/
static int mes_process_win_enabled (void)
{
  mes_t *mpc = mes_process_get ();
  if (mpc)
    return mpc->win_enabled;
  return (1);
}                               /* mes_process_win_enabled */

/*---------------------------------------------------------------------------*/
char *mes_get_std_path (void)
{
  mes_t *mpc = mes_process_get ();
  if (mpc)
    return mpc->std_path;
  return (NULL);
}                               /* mes_get_std_path */

/*---------------------------------------------------------------------------*/
static void mes_process_alloc ()
{
  mes_t *mpi;
  int tid = getthreadid ();
  if (mes_process_get ())
    return;
  if (mes_process_n >= MES_PROCESS_MAX)
    return;
  if (!mes_process)
    return;
  else {
    int i, j = -1;
    for (i = mes_process_n - 1; i >= 0; i--)
      if (!mes_process[i]) {
        j = i;
        break;
      }
    if (i == -1 && mes_process_n < MES_PROCESS_MAX) {
      mpi = (mes_t *) calloc (1, sizeof (mes_t));
      if (!mpi)
        return;
      mes_process[mes_process_n++] = mpi;
    }
    else if (j >= 0) {
      mpi = (mes_t *) calloc (1, sizeof (mes_t));
      if (!mpi)
        return;
      mes_process[j] = mpi;
    }
    else
      return;
  }

  mpi->thread_id = tid;
  mpi->std_path = NULL;
  mpi->logfile = NULL;
  mpi->argc = 0;
  mpi->argv = NULL;

  mpi->win_enabled = 1;
  mpi->enabled = 1;
  mpi->win_fkt = NULL;
}                               /* mes_process_alloc */

/*---------------------------------------------------------------------------*/
static void mes_aux_va (int flags, char *format, va_list args)
{
  char str[1024];
  char *tmp;

  tmp = mprintf_va_dyn (str, sizeof (str), format, args);
  if (!tmp)
    return;
  mes_smart (flags, tmp, -1);
  if (tmp - str)
    free (tmp);
  return;
}                               /* mes_aux_va */

/*---------------------------------------------------------------------------*/
static void mes_aux (int flags, char *format, ...)
{
  va_list args;
  va_start (args, format);
  mes_aux_va (flags, format, args);
  return;
}                               /* mes_aux */

/*---------------------------------------------------------------------------*/
static int mes_filename_check (const char *filename)
{
  int len;

  if (!filename)
    return (-1);
  len = strlen (filename);
  if (len <= 0)
    return (-1);
#if defined(_WINDOWS)
  /* WORK AROUND !! 
     The following code line is included to catch an error of the 
     network driver of Windows NT. 
   */

  if (filename[len - 1] == '\\')
    return (-1);

  /* End of the work around */
#endif
  return (0);
}                               /* mes_filename_check */

/*---------------------------------------------------------------------------*/
static void mes_arg_free (void)
{
  mes_t *mpc = mes_process_get ();
  if (!mpc)
    return;
  if (mpc->argv) {
    while (mpc->argc > 0) {
      mpc->argc--;
      if (mpc->argv[mpc->argc])
        free (mpc->argv[mpc->argc]);
    }
    free (mpc->argv);
    mpc->argv = NULL;
  }
  mpc->argc = 0;
}                               /* mes_arg_free */

/*----------------------------------------------------------------------------*/
static void mes_init_std_path (char *logfile)
{
  mes_t *mpc = mes_process_get ();
  int i;

  if (mpc) {
    if (mpc->std_path)
      free (mpc->std_path);
    for (i = logfile ? strlen (logfile) : 0;
         i && logfile[i - 1] - PathSepChar; i--);
    mpc->std_path = malloc (i + 1);
    if (mpc->std_path) {
      memcpy (mpc->std_path, logfile, i);
      mpc->std_path[i] = '\0';
    }
  }
}                               /* mes_init_std_path */



/*============================================================================*/
void mes_smart (int flags, const char *txt, int bytes)
{
  char tmp[2] = { 0 };
  int slen;

  if (!mes_process_enabled ())
    return;
  if (flags & MES_FLAG_TIME)
    mes_time ();
  if (!txt)
    return;
  if (bytes < 0) {
    slen = strlen (txt);
    bytes = slen;
  }
  else {
    char *p = memchr (txt, 0, bytes);
    if (p)
      slen = txt - p;
    else
      slen = bytes + 1;
  }
  if (bytes <= 0)
    return;
  if (slen > bytes) {
    tmp[0] = txt[bytes - 1];
    tmp[bytes - 1] = 0;
  }

  if (flags & (MES_FLAG_FILE | MES_FLAG_TIME)) {
    FILE *fmes;
    int size;
    char *logfile = NULL;
    mes_t *mpc = mes_process_get ();
    if (mpc)
      logfile = mpc->logfile;

    if (!logfile)
      logfile = mes_default_logfile;
    fmes = fopen (logfile, "rb");
    if (fmes) {
      fseek (fmes, 0, SEEK_END);
      size = ftell (fmes);
      fclose (fmes);
      if (size > MAX_MES_FILE_LENGTH) {
        char bakfile[300];
        if (mpc && mpc->std_path)
          sprintf (bakfile, "%s%s", mpc->std_path, mes_default_bak_logfile);
        else
          sprintf (bakfile, "%s", mes_default_bak_logfile);
        remove (bakfile);
        if (rename (logfile, bakfile)) {
          mes_aux (MES_FLAG_WIN,
                   "\nFehler: Kann Datei %s nicht in %s umbenennen.\n\n",
                   logfile, bakfile);
        }
      }
    }

    fmes = fopen (logfile, "at");
    if (fmes) {
      fputs (txt, fmes);
      fputs (tmp, fmes);
      fclose (fmes);
    }
  }
  if ((flags & MES_FLAG_WIN) && mes_process_win_enabled ()) {
    mes_t *mpc = mes_process_get ();
    if (mpc && mpc->win_fkt) {
      (*mpc->win_fkt) (txt);
      (*mpc->win_fkt) (tmp);
    }

    else {
      fputs (txt, stdout);
      fputs (tmp, stdout);
      fflush (stdout);
    }
  }
  if (slen > bytes)
    ((char *) (txt))[bytes - 1] = tmp[0];
}                               /* mes_smart */

/*============================================================================*/
void mes_exit (void)
{
  int i, tid = getthreadid ();
  for (i = mes_process_n - 1; i >= 0; i--) {
    mes_t *mpi = mes_process[i];
    if (mpi)
      if (mpi->thread_id == tid) {
        if (mpi->logfile)
          free (mpi->logfile);
        if (mpi->std_path)
          free (mpi->std_path);
        mes_arg_free ();
        free (mpi);
        mes_process[i] = NULL;
        return;
      }
  }
}                               /* mes_exit */

/*============================================================================*/
void mes_init_args (int argc, char *argv[])
{
  mes_t *mpc = mes_process_get ();
  if (!argv || argc < 1)
    return;
  mes_arg_free ();
  mpc->argv = calloc (1, argc * sizeof (*mpc->argv));
  if (mpc->argv)
    while (mpc->argc < argc && *argv) {
      int len = *argv ? strlen (*argv) + 1 : 1;
      mpc->argv[mpc->argc] = malloc (len);
      if (!mpc->argv[mpc->argc])
        break;
      if (*argv)
        memcpy (mpc->argv[mpc->argc], *argv, len);
      else
        mpc->argv[mpc->argc][0] = 0;
      argv++;
      mpc->argc++;
    }
}                               /* mes_init_args */

/*============================================================================*/
void mes_init_logfile (char *logfile)
{
  FILE *tst;
  mes_t *mpc = mes_process_get ();

  if (!logfile)
    return;

  mes_init_std_path (logfile);

  tst = fopen (logfile, "at");
  if (tst) {
    fclose (tst);
    if (mpc) {
      int len = strlen (logfile);
      if (mpc->logfile)
        free (mpc->logfile);
      mpc->logfile = (char *) malloc (len + 1);
      if (mpc->logfile) {
        strcpy (mpc->logfile, logfile);
        mpc->logfile[len] = '\0';
      }
      return;
    }
  }
}                               /* mes_init_logfile */

/*============================================================================*/
void mes_init_winfct (void (*winfct) (const char *))
{
  mes_t *mpc = mes_process_get ();
  if (mpc && winfct)
    mpc->win_fkt = winfct;
}                               /* mes_init_winfct */

/*============================================================================*/
void mes_init (char *logfile, void (*winfct) (const char *), int argc,
               char *argv[])
{
  mes_process_alloc ();
  mes_init_args (argc, argv);
  mes_init_logfile (logfile);
  mes_init_winfct (winfct);
}                               /* mes_init */

/*============================================================================*/
int mes_win_ability (int on)
{
  mes_t *mpc = mes_process_get ();
  if (mpc) {
    int res = mpc->win_enabled;
    mpc->win_enabled = !!on;
    return (res);
  }
  return (1);
}                               /* mes_win_ability */

/*============================================================================*/
int mes_ability (int on)
{
  mes_t *mpc = mes_process_get ();
  if (mpc) {
    int res = mpc->enabled;
    mpc->enabled = !!on;
    return (res);
  }
  return (1);
}                               /* mes_ability */

/*============================================================================*/
void mes_time (void)
{
  mes_t *mpc = mes_process_get ();
  time_t now = time (NULL);
  char *timestr = ctime (&now);
  int len = strlen (timestr);
  char txt[256];

  if (!len)
    return;
  if (timestr[len - 1] == '\n')
    timestr[len - 1] = ' ';
  mes_file ("\n***** ");
  if (mpc) {
    if (mpc->argc == 1) {
      mes_file (mpc->argv[0]);
      mes_file (":");
    }
    else if (mpc->argc > 1) {
      int i;
      mes_file (mpc->argv[0]);
      mes_file ("(");
      mes_file (mpc->argv[1]);
      for (i = 2; i < mpc->argc; i++) {
        mes_file (",");
        mes_file (mpc->argv[i]);
      }
      mes_file (")");
      mes_file (":");
    }
  }
  mes_file (timestr);
  sprintf (txt, "(%.2f sec)", clock () / (float) CLOCKS_PER_SEC);
  mes_file (txt);
  mes_file (" *****:\n");
}                               /* mes_time */


/*============================================================================*/
void mes_va (int flags, int line, char *xproc, char *proc, char *format,
             va_list args)
{
  char digstr[32] = { 0 };

  if (!format && !xproc && !proc)
    return;
  if (line + 1)
    mprintf (digstr, sizeof (digstr), "(%u)", line);
  if (flags & MES_FLAG_TIME) {
    mes_time ();
    flags &= ~MES_FLAG_TIME;
    flags |= MES_FLAG_FILE;
  }
  if (!xproc)
    xproc = proc;
  if (!proc)
    proc = xproc;
  if (proc) {
    if (flags & MES_FLAG_FILE)
      mes_file (xproc);
    if (flags & MES_FLAG_WIN)
      mes_win (proc);
    mes_smart (flags, digstr, -1);
    if (format)
      mes_smart (flags, ": ", -1);
  }

  if (!format)
    mes_smart (flags, "\n", -1);
  else
    mes_aux_va (flags, format, args);
}                               /* mes_va */

/*============================================================================*/
void mes (int flags, int line, char *xproc, char *proc, char *format, ...)
{
  va_list args;
  va_start (args, format);
  mes_va (flags, line, xproc, proc, format, args);
}                               /* mes */

/*============================================================================*/
void mes_printf (int prot_flags, char *format, ...)
{
  va_list args;
  char *txt;

  if (!prot_flags)
    return;
  if (!format) {
    mes_time ();
    mes_file ("Call of mes_printf without format string");
    return;
  }
  va_start (args, format);
  txt = mprintf_va (NULL, 0, format, args);
  if (!txt) {
    mes_time ();
    mes_file_win ("Call of mes_printf with format string\"");
    mes_file_win (format);
    mes_file_win ("\" without success\n");
    return;
  }
  if (prot_flags & MES_FLAG_TIME)
    mes_time ();
  mes_smart (prot_flags, txt, -1);
  free (txt);
}                               /* mes_printf */

/*============================================================================*/
void mes_fformat (char *txt, char *logfile, int line, char *proc_info)
{
  mes_time ();
  if (proc_info && strlen (proc_info)) {
    mes_smart (MES_FLAG_FILE, proc_info, -1);
    mes_file (":");
  }
  mes_file_win ("format error");
  if (logfile && strlen (logfile)) {
    mes_file_win (" in file ");
    mes_file_win (logfile);
  }
  if (line >= 0)
    mes_aux (MES_FLAG_FILE_WIN, ": line %d", line);
  if (logfile && strlen (txt)) {
    mes_file_win (" (");
    mes_file_win (txt);
    mes_file_win (")\n");
  }
  else
    mes_file_win ("\n");

}                               /* mes_fformat */


/*============================================================================*/
void mes_err (char *txt, int error_nr, char *proc_info)
{
  mes_time ();
  if (proc_info && strlen (proc_info)) {
    mes_file_win (proc_info);
    mes_file_win (":");
  }
  if (error_nr >= 0
      && error_nr < sizeof (mes_err_txt) / sizeof (*mes_err_txt)) {
    mes_file_win (mes_err_txt[error_nr]);
  }
  if (txt) {
    mes_file_win (" (");
    mes_file_win (txt);
    mes_file_win (")\n");
  }
  else
    mes_file_win ("\n");

}                               /* mes_err */

/*============================================================================*/
void mes_proc_start (char *proc_info)
{
  mes_t *mpc = mes_process_get ();
  int i;

  mes_time ();
  if (proc_info) {
    mes_file (proc_info);
    mes_file (":");
  }
  mes_file (" ***** PROGRAM STARTED *****\n");
  if (mpc)
    for (i = 0; i < mpc->argc; i++) {
      if (!i)
        mes_file ("program call name is : ");
      else
        mes_aux (MES_FLAG_FILE, "parameter %10d : ", i);
      mes_file (mpc->argv[i]);
      mes_file ("\n");
    }

}                               /* mes_proc_start */

/*============================================================================*/
int mes_insert (FILE * fp, char src, int cnt)
{
  int i = cnt;

  if (fp && fp - stdout)
    for (i = 0; i < cnt && !mes_fputc (fp, src); i++);
  else
    for (i = 0; i < cnt; i++)
      mes_smart (MES_FLAG_WIN, &src, 1);
  if (i == cnt)
    return (0);
  return (-1);
}                               /* mes_insert */


/******************************************************************************/
/******************************************************************************/

/*============================================================================*/
void *mes_malloc (int bytes)
{
  void *res;

  if (bytes <= 0)
    bytes = 1;
  res = malloc (bytes);
  if (res)
    return (res);
  else
    mes_aux (MES_FLAG_TIME_WIN, "malloc: could not allocate %d bytes\n",
             bytes);
  return (NULL);
}                               /* mes_malloc */

/*============================================================================*/
void *mes_calloc (int bytes)
{
  void *res;

  if (bytes <= 0)
    bytes = 1;
  res = calloc (1, bytes);
  if (res)
    return (res);
  else
    mes_aux (MES_FLAG_TIME_WIN, "calloc: could not allocate %d bytes\n",
             bytes);
  return (NULL);
}                               /* mes_calloc */

/*============================================================================*/
int mes_realloc (void **mem, int bytes)
{
  void *res;

  if (bytes <= 0)
    bytes = 1;
  if (!mem)
    return (-1);
  if (!*mem)
    res = malloc (bytes);
  else
    res = realloc (*mem, bytes);
  if (res) {
    *mem = res;
    return (0);
  }
  else
    mes_aux (MES_FLAG_TIME_WIN,
             "realloc: could not reallocate %d bytes\n", bytes);
  return (-1);
}                               /* mes_realloc */

/*============================================================================*/
FILE *mes_fopen (const char *filename, char *attrstr)
{
  FILE *fp;

  if (mes_filename_check (filename))
    goto STOP;
  if (!attrstr)
    goto STOP;
  if (!strcmp (filename, "stdout")) {
    return (stdout);
  }
  fp = fopen (filename, attrstr);
  if (fp) {
    if (!strchr (attrstr, 'b') && !strchr (attrstr, 't')) {
      mes_file_win ("fopen: file \"");
      mes_file_win (filename);
      mes_file_win ("\" opened with ambiguous attributes \"");
      mes_file_win (attrstr);
      mes_file_win ("\"\n");
    }
    return (fp);
  }

STOP:
  mes_time ();
  mes_file_win ("fopen: could not open file \"");
  mes_file_win (filename);
  mes_file_win ("\" with attribute \"");
  mes_file_win (attrstr);
  mes_file_win ("\"\n");
  return (NULL);
}                               /* mes_fopen */

/*============================================================================*/
int mes_fread_quiet (FILE * fp, void *mem, int bytes)
{
  if (!bytes)
    return (0);
  if (mem && fp)
    return (fread (mem, 1, bytes, fp));
  else
    mes_aux (MES_FLAG_TIME_WIN,
             "fread: could not read %d bytes from FILE(%p) to mem(%p)\n",
             bytes, fp, mem);
  return (-1);
}                               /* mes_fread_quiet */

/*============================================================================*/
int mes_fread (FILE * fp, void *mem, int bytes)
{
  if (!bytes)
    return (0);
  if (mem && fp && fread (mem, 1, bytes, fp) == (unsigned int) bytes)
    return (0);

  else
    mes_aux (MES_FLAG_TIME_WIN,
             "fread: could not read %d bytes from FILE(%p) to mem(%p)\n",
             bytes, fp, mem);
  return (-1);
}                               /* mes_fread */

/*============================================================================*/
int mes_fwrite (FILE * fp, void *mem, int bytes)
{
  if (!fp || !mem)
    bytes = -1;
  else if (bytes < 0)
    bytes = strlen (mem);
  if (!bytes)
    return (0);
  if (bytes > 0 && fwrite (mem, 1, bytes, fp) == (unsigned int) bytes)
    return (0);
  else
    mes_aux (MES_FLAG_TIME_WIN,
             "fwrite: could not write %d bytes from mem(%p) to FILE(%p)\n",
             bytes, mem, fp);
  return (-1);
}                               /* mes_fwrite */


/*============================================================================*/
int mes_fputc (FILE * fp, char chr)
{
  if (fp && fputc (chr, fp) != EOF)
    return (0);

  else
    mes_aux (MES_FLAG_TIME_WIN,
             "fputc: could not write byte %X to FILE(%p)\n",
             (int) chr & 0xFF, fp);
  return (-1);
}                               /* mes_fputc */

/*============================================================================*/
int mes_fputs (FILE * fp, char *str)
{
  if (fp && str && fputs (str, fp) != EOF)
    return (0);

  if (str)
    mes_aux (MES_FLAG_TIME_WIN, "fputs: could not write string %s\n", str);
  else
    mes_aux (MES_FLAG_TIME_WIN, "fputs: could not write 0 pointer\n");

  return (-1);
}                               /* mes_fputs */

/*============================================================================*/
int mes_fgetc (FILE * fp)
{
  int res = fp ? fgetc (fp) : EOF;
  if (res != EOF)
    return (res);
  else
    mes_aux (MES_FLAG_TIME_WIN, "fgetc: end of FILE(%p)\n", fp);
  return (res);
}                               /* mes_fgetc */


/*============================================================================*/
int mes_fflush (FILE * fp)
{
  int res = fp ? fflush (fp) : -1;

  if (!res)
    return (0);
  else
    mes_aux (MES_FLAG_TIME_WIN, "fflush: could not flush FILE(%p)\n", fp);
  return (res);
}                               /* mes_fflush */

/*============================================================================*/
int mes_fprintf (FILE * fp, char *format, ...)
{
  va_list args;
  char *txt;

  if (!format)
    return (0);

  va_start (args, format);
  txt = mprintf_va (NULL, 0, format, args);
  if (!txt) {
    mes_time ();
    mes_file_win ("sprintf_va: call with format string\"");
    mes_file_win (format);
    mes_file_win ("\" without success\n");
    return (-1);
  }
  if (fp && fp - stdout)
    mes_fputs (fp, txt);
  else
    mes_win (txt);
  free (txt);
  return (1);
}                               /* mes_fprintf */

/*============================================================================*/
int mes_fseek (FILE * fp, long offset, int fromwhere)
{
  int res = fp ? fseek (fp, offset, fromwhere) : -1;

  if (!res)
    return (0);
  mes_aux (MES_FLAG_TIME_WIN, "fseek: could not position FILE(%p) at %ld",
           fp, offset);
  switch (fromwhere) {
  case SEEK_SET:
    mes_aux (MES_FLAG_FILE_WIN, "\n");
    break;
  case SEEK_END:
    mes_aux (MES_FLAG_FILE_WIN, " from the end\n");
    break;
  case SEEK_CUR:
    mes_aux (MES_FLAG_FILE_WIN, " from current position\n");
    break;
  default:
    mes_aux (MES_FLAG_FILE_WIN, " with undefinded offset %d\n", fromwhere);
    break;
  }
  return (res);
}                               /* mes_fseek */

/*============================================================================*/
#ifdef WIN32
int mes_fseek64 (FILE * fp, unsigned int uoff, unsigned int loff,
                 int fromwhere)
{
/* lower 32bit in lpos, upper 32bit in upos */
  int fh;                       /* file handle */
  __int64 pos64;
  __int64 off64;
  if (!fp)
    return (0);

  fh = _fileno (fp);
  off64 = (((__int64) uoff) << 32) + (__int64) loff;
  pos64 = _lseeki64 (fh, off64, fromwhere);
  if (pos64 == -1L) {
    mes_file_win ("_lseeki64 failed:");
    mes_file_win (strerror (errno));
    return (0);
  }
  return (1);
}

/*============================================================================*/
int mes_ftell64 (int fh, unsigned int *upos, unsigned int *lpos)
{
  /* lower 32bit in *lpos, upper 32bit in *upos */
  __int64 pos64;
  if (fh < 0)
    return (0);

  pos64 = _lseeki64 (fh, 0, SEEK_CUR);
  if (pos64 == -1L) {
    mes_file_win ("_lseeki64 failed:");
    mes_file_win (strerror (errno));
    return (0);
  }
  *upos = (unsigned int) (pos64 >> 32);
  *lpos = (unsigned int) pos64;
  return (1);
}
#endif /* WIN32 */

/*============================================================================*/
int mes_ftell (FILE * fp)
{
  int res = fp ? ftell (fp) : -1;

  if (res != -1)
    return (res);
  else
    mes_aux (MES_FLAG_TIME_WIN,
             "ftell: could not find current position of FILE(%p)\n", fp);
  return (res);
}                               /* mes_ftell */

/*============================================================================*/
int mes_remove (char *filename)
{
#define CUR_PROC "mes_remove"
  int res = -1;

  if (mes_filename_check (filename))
    goto STOP;
  res = remove (filename);
  if (!res)
    return (0);
STOP:
  mes_time ();
  mes_file_win ("remove: could not remove file \"");
  mes_file_win (filename);
  mes_file_win ("\";");
  if (res != -1)
    mes_file_win (strerror (errno));
  mes_file_win ("\n");
  return (res);
#undef CUR_PROC
}                               /* mes_remove */


/*============================================================================*/
int mes_rename (char *oldname, char *newname)
{
#define CUR_PROC "mes_rename"
  int res = -1;

  if (mes_filename_check (oldname))
    goto STOP;
  if (mes_filename_check (newname))
    goto STOP;

#if defined(_PPC_) || defined(WIN32)
  if (strcmp (oldname, newname)) {
    FILE *fp;
    fp = fopen (newname, "rb");
    if (fp) {
      fclose (fp);
      mes_remove (newname);
    }
  }
  else
    return (0);
#endif
  res = rename (oldname, newname);
  if (!res)
    return (res);

STOP:
  mes_time ();
  mes_file_win ("rename: could not rename \"");
  mes_file_win (oldname);
  mes_file_win ("\" -> \"");
  mes_file_win (newname);
  mes_file_win ("\"; ");
  if (res != -1)
    mes_file_win (strerror (errno));
  mes_file_win ("\n");
  return (res);
#undef CUR_PROC
}                               /* mes_rename */

/*============================================================================*/
int mes_move (char *oldname, char *newname)
{
#define CUR_PROC "mes_move"
  int res = -1;
  int tmp;

  if (mes_filename_check (oldname))
    goto STOP;
  if (mes_filename_check (newname))
    goto STOP;
  if (!strcmp (oldname, newname))
    goto STOP;

  /* first try to rename */
  tmp = mes_ability (0);
  if (!mes_rename (oldname, newname)) {
    res = 0;
    mes_ability (tmp);
    goto STOP;
  }
  mes_ability (tmp);

  /* need to copy */
  if (mes_copy (oldname, newname))
    goto STOP;
  mes_remove (oldname);

  res = 0;
STOP:
  if (res < 0) {
    mes_time ();
    mes_file_win ("move: could not move ");
    mes_file_win (oldname);
    mes_file_win (" -> ");
    mes_file_win (newname);
    mes_file_win ("\n");
  }
  return (res);
#undef CUR_PROC
}                               /* mes_move */


/*============================================================================*/
int mes_copy (char *oldname, char *newname)
{
#define CUR_PROC "mes_copy"
  int res = -1;
  FILE *dst = NULL;
  FILE *src = NULL;
  char *buf = NULL;

  if ((dst = mes_fopen (newname, "wb")) == NULL)
    goto STOP;
  if ((src = mes_fopen (oldname, "rb")) == NULL)
    goto STOP;
  if ((buf = mes_malloc (0x10000 * sizeof (char))) == NULL)
    goto STOP;

  while (!feof (src)) {
    int cnt = fread (buf, 1, 0x10000, src);
    if (cnt <= 0)
      break;
    if (fwrite (buf, 1, cnt, dst) <= 0)
      goto STOP;
  }
  res = 0;
STOP:
  if (buf)
    free (buf);
  if (src)
    fclose (src);
  if (dst)
    fclose (dst);
  return (res);
# undef CUR_PROC
}                               /* mes_copy */

/*============================================================================*/
FILE *mes_tmpfopen (char *path)
{
# define CUR_PROC "mes_tmpfopen"
  FILE *fp;
  char name[16];
  char tmpname[L_tmpnam + 16];
  int i;

  if (path)
    strncpy (tmpname, path, L_tmpnam);
  else
    tmpname[0] = 0;

  for (i = 0; i < 0x10000; i++) {
    sprintf (name, "%80X.TMP", rand ());
    strcat (tmpname, name);
    fp = fopen (tmpname, "rb");
    if (fp) {
      fclose (fp);
      continue;
    }
    fp = fopen (tmpname, "w+b");
    if (fp)
      return (fp);
    break;
  }
  mes_time ();
  mes_file_win ("tmpfopen: no success\n");
  return (NULL);
# undef CUR_PROC
}                               /* mes_tmpfopen */


/*============================================================================*/
int mes_tmpfclose (FILE ** fp)
{
  if (!fp)
    return (0);
  if (!*fp)
    return (0);
  fclose (*fp);
  (*fp) = 0;
  return (0);
}                               /* mes_tmpfclose */


/*============================================================================*/
FILE *mes_tmpfile (void)
{
# define CUR_PROC "mes_tmpfile"
  FILE *fp = tmpfile ();
  if (fp)
    return (fp);
  else {
    mes_time ();
    mes_file_win ("tmpfile: no success\n");
  }
  return (NULL);
# undef CUR_PROC
}                               /* mes_tmpfile */
