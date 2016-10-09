#!gmake

VER=1.6
NM=GemSIM_v$(VER)
AR=$(NM).tar.gz
URL=http://downloads.sourceforge.net/project/gemsim/$(AR)

$(NM): $(AR)
	tar zxvf $<

$(AR):
	wget $(URL)

.PHONY: clean
clean:
	rm -fr $(NM)
	rm -f $(AR)
