#!gmake

VER=1.5.1
NM=art_Illumina_src
NM_VER=$(NM)-$(VER)
AR=$(NM_VER).tar.gz
URL=http://www.niehs.nih.gov/research/resources/assets/docs/$(AR)

#$(NM)/???: $(NM)
#	cd $(NM) && perl Makefile.PL --defaultdeps
#	cd $(NM) && make

$(NM): $(AR)
	tar zxvf $<

$(AR):
	wget $(URL)
