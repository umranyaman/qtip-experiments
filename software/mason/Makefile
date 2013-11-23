#
# Mason simulator - assumues cmake & gcc tools are already installed
#

SVN=$(shell which svn)
MASON_SRC=seqan-trunk/build/Release/bin/mason
SITE=http://svn.seqan.de

mason: $(MASON_SRC)
	ln -s -f seqan-trunk/build/Release/bin/mason

$(MASON_SRC): seqan-trunk/build
	cd seqan-trunk/build && mkdir -p Release && cd Release && cmake -DCMAKE_BUILD_TYPE=Release ../.. && $(MAKE) mason

seqan-trunk/build: $(SVN)
	$(SVN) co $(SITE)/seqan/trunk seqan-trunk
	mkdir -p $@
