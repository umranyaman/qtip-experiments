#!gmake

# Check for new versions here:
# https://github.com/bmvdgeijn/WASP/releases

NM=WASP
VER=0.2.1
AR=v$(VER).tar.gz
FULL_NM=$(NM)-$(VER)
URL=https://github.com/bmvdgeijn/$(NM)/archive/$(AR)

$(FULL_NM)/mapping/Snakefile: $(AR)
	tar zvxf $<
	touch $@

$(AR):
	wget $(URL)
