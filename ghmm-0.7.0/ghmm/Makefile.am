#*******************************************************************************
#
#       This file is part of the General Hidden Markov Model Library,
#       GHMM version __VERSION__, see http://ghmm.org
#
#       Filename: ghmm/ghmm/Makefile.am
#       Authors:  Wasinee Rungsarityotin, Benjamin Georgi, Alexander Schliep, Benjamin Rich
#
#       Copyright (C) 1998-2004 Alexander Schliep
#       Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
#       Copyright (C) 2002-2004 Max-Planck-Institut fuer Molekulare Genetik,
#                               Berlin
#
#       Contact: schliep@ghmm.org
#
#       This library is free software; you can redistribute it and/or
#       modify it under the terms of the GNU Library General Public
#       License as published by the Free Software Foundation; either
#       version 2 of the License, or (at your option) any later version.
#
#       This library is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#       Library General Public License for more details.
#
#       You should have received a copy of the GNU Library General Public
#       License along with this library; if not, write to the Free
#       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
#
#       This file is version $Revision: 1269 $
#                       from $Date: 2005-08-18 17:10:19 +0200 (Thu, 18 Aug 2005) $
#             last change by $Author: grunau $.
#
#******************************************************************************/

SUBDIRS = .
MAINTAINERCLEANFILES = Makefile.in

INCLUDES = -I$(top_srcdir)
lib_LTLIBRARIES = libghmm.la

libghmm_la_LDFLAGS = -version-info @GHMM_VERSION_CURRENT@:@GHMM_VERSION_REVISION@:@GHMM_VERSION_AGE@
libghmm_la_SOURCES= foba.c foba.h \
                  matrix.c  matrix.h \
                  model.c model.h\
                  modelutil.c modelutil.h \
                  cluster.c cluster.h \
                  scluster.c scluster.h \
                  root_finder.c root_finder.h \
                  reestimate.c reestimate.h \
                  randvar.c randvar.h \
                  rng.c rng.h \
                  sequence.c sequence.h \
                  vector.c vector.h \
                  viterbi.c viterbi.h \
                  sfoba.c sfoba.h \
                  sviterbi.c sviterbi.h \
                  smodel.c smodel.h \
                  sdmodel.c sdmodel.h \
                  sdviterbi.c \
                  sdfoba.c sdfoba.h \
                  ghmm.h \
                  sreestimate.c sreestimate.h \
                  sgenerate.c sgenerate.h \
                  smixturehmm.c smixturehmm.h \
                  const.h \
                  mes.c mes.h \
                  mprintf.c mprintf.h \
                  scanner.c scanner.h \
                  gauss_tail.c gauss_tail.h \
                  smap_classify.c smap_classify.h\
                  kbestbasics.c kbestbasics.h \
                  kbest.c kbest.h \
                  gradescent.c gradescent.h \
                  discrime.c discrime.h \
		  psequence.c psequence.h \
		  pmodel.c pmodel.h \
		  linkedlist.c linkedlist.h \
		  pviterbi.c pviterbi.h \
		  pviterbi_propagate.c pviterbi_propagate.h


pkginclude_HEADERS = ghmm.h \
		  foba.h \
                  model.h \
                  cluster.h \
                  scluster.h \
                  root_finder.h \
                  reestimate.h \
                  sequence.h \
                  viterbi.h \
                  sfoba.h \
                  sviterbi.h \
                  smodel.h \
		  sdmodel.h \
                  sdfoba.h \
                  sreestimate.h \
                  sgenerate.h \
                  smixturehmm.h \
                  scanner.h \
                  gauss_tail.h \
                  smap_classify.h \
                  kbest.h \
                  gradescent.h \
                  discrime.h \
		  psequence.h \
		  pmodel.h \
		  linkedlist.h \
		  pviterbi.h \
		  pviterbi_propagate.h \
		  vector.h \
		  matrix.h \
		  const.h \
		  rng.h \
		  mes.h \
		  randvar.h \
		  internal.h


# must be eleminated...
pkgdata_DATA = PHI_001_20.dat
EXTRA_DIST= PHI_001_20.dat

COPYRIGHT = $(top_srcdir)/COPYRIGHT.txt

dist-hook:
	test -z "$(AWK)" && exit 0 ;\
	TMPFILE="$(distdir)/tmp"; \
	for f in $(DISTFILES) ; do \
	  distfile=$(distdir)/$$f ; \
	  if grep '__copyright__ *$$' $$distfile >/dev/null ; then \
	    mv $$distfile $$TMPFILE; \
	    $(AWK) 'BEGIN {nr=0;while(getline < "$(COPYRIGHT)" ==1) copyright[nr++]=$$0;} /__copyright__ *$$/{sub("__copyright__.*","");for(i=0;i<nr;i++) print $$0 copyright[i];} //' \
	    $$TMPFILE >$$distfile; \
	    touch -m -a -r $$TMPFILE $$distfile; \
	  fi; \
	done ; \
	rm -f $$TMPFILE

