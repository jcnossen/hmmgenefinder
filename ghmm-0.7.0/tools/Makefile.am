#  author       : Achim Gaedke
#  filename     : /zpr/bspk/src/hmm/ghmm/tools/Makefile.am
#  created      : DATE: April 2001
#  $Id: Makefile.am 221 2002-03-05 10:03:49Z pipenb $
#
# Copyright (C) 1998-2005 Alexander Schliep
# Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
# Copyright (C) 2002-2005 Max-Planck-Institut fuer Molekulare Genetik, Berlin
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
# 
# 
# 

SUBDIRS = .
DISTCLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

BUILT_SOURCES = 
INCLUDES = -I$(top_srcdir)
bin_PROGRAMS = probdist cluster scluster smix_hmm

probdist_SORUCES = probdist.c
cluster_SOURCES  = cluster.c
scluster_SOURCES = scluster.c
smix_hmm_SOURCES = smix_hmm.c

LDADD = $(top_builddir)/ghmm/.libs/libghmm.a
bin_SCRIPTS = ghmm-config

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








