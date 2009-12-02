#  author       : Achim Gaedke
#  filename     : ghmm/tests/Makefile.am
#  created      : DATE: April 2001
#  $Id: Makefile.am 1214 2005-07-06 12:14:32Z grunau $
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

CCLD     = $(CC)   -lstdc++

SUBDIRS = data

# tests and examples

AUTOMAKE_OPTIONS = dejagnu
EXPECT=$(EXPECT_FOR_DEJAGNU)
RUNTEST=$(RUNTEST_FOR_DEJAGNU)

CLEANFILES = ghmm.log ghmm.sum site.exp site.bak
MAINTAINERCLEANFILES = Makefile.in
# stuff from dejagnu

INCLUDES=-I$(top_srcdir)

#these tests will not be installed
# in progress: test_sdfoba nullmodel ciscreen test_sdmodel 
noinst_PROGRAMS = root_finder_test \
                  coin_toss_test \
                  two_states_three_symbols \
                  chmm \
                  shmm_viterbi_test \
                  generate_PHI \
                  test_gsl_ran_gaussian_tail \
                  sequences_test \
                  sequences_old_format \
                  label_higher_order_test

LDADD = $(top_builddir)/ghmm/.libs/libghmm.a

COPYRIGHT = $(top_srcdir)/COPYRIGHT.txt

# creates copyright tags and copies some files
dist-hook:
	DEJAGNU_FILES=`(cd $(srcdir) && find ghmm.tests lib config -name \*.exp)` ;\
	for f in $$DEJAGNU_FILES ; do \
	  dirpart=`dirname $$f`; \
	  test -d $(distdir)/$$dirpart || mkdir $(distdir)/$$dirpart ; \
	  cp -p $(srcdir)/$$f $(distdir)/$$dirpart ; \
	done; \
	test -z "$(AWK)" && exit 0 ;\
	TMPFILE="$(distdir)/tmp" ; \
	for f in $(DISTFILES) $$DEJAGNU_FILES ; do \
	  distfile=$(distdir)/$$f ; \
	  if grep '__copyright__ *$$' $$distfile >/dev/null ; then \
	    mv $$distfile $$TMPFILE; \
	    $(AWK) 'BEGIN {nr=0;while(getline < "$(COPYRIGHT)" ==1) copyright[nr++]=$$0;} /__copyright__ *$$/{sub("__copyright__.*","");for(i=0;i<nr;i++) print $$0 copyright[i];} //' \
	    $$TMPFILE >$$distfile; \
	    touch -m -a -r $$TMPFILE $$distfile; \
	  fi; \
	done ; \
	rm -f $$TMPFILE
