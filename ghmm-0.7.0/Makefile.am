#  author       : Achim Gaedke
#  filename     : ghmm/Makefile.am
#  created      : DATE: April 2001
#  $Id: Makefile.am 360 2003-12-19 19:21:14Z schliep $
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

SUBDIRS = @SUBDIRS@

MAINTAINERCLEANFILES = \
aclocal.m4 \
configure \
Makefile.in \
stamp-h.in \
config.h.in \
mkinstalldirs \
missing \
install-sh \
INSTALL \
COPYING \
$(distdir).tar.gz -r $(distdir) \
ltmain.sh \
config.guess \
config.sub \
$(distdir).tar.gz \
-r $(distdir)

EXTRA_DIST = \
COPYRIGHT.txt \
GPL.txt \
LGPL.txt \
ToDo \
win_config.h \
autogen.sh \
README.Win

COPYRIGHT = $(top_srcdir)/COPYRIGHT.txt

dist-hook:
	test -z "$(AWK)" && exit 0 ;\
	TMPFILE="$(distdir)/tmp" ; \
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
