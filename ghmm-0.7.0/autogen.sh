#! /bin/sh
#  author       : Achim Gaedke
#  filename     : ghmm/autogen.sh
#  created      : DATE: April 2001
#  $Id: autogen.sh 1135 2005-05-11 14:42:45Z schliep $

#move GNU m4 to head of PATH
sep_path=`echo $PATH|sed -n 's/:/ /gp'`
for f in $sep_path ; do
  test -f $f/m4 && echo `eval $f/m4 --version 2>/dev/null` |grep GNU >/dev/null && PATH=$f:$PATH
done

#creates ltmain.sh (Call first!!!)
libtoolize --automake --force --copy

#makes aclocal.m4 from acinclude.m4 and other files
aclocal

#scans configure.in and creates config.h.in
autoheader 

#creates Makefile.in from Makefile.am
automake --add-missing --force-missing

#creates configure from configure.in
autoconf


