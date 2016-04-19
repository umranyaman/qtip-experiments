#!gmake

VER=1.0beta.18
AR=v$(VER).tar.gz
URL=https://github.com/amplab/snap/archive/$(AR)
NM=snap
NM_VER=snap-$(VER)

snap-aligner: $(NM)/snap-aligner
	ln -s -f $< $@

$(NM)/snap-aligner: $(NM)/Makefile
	make -C $(NM) # CXX=g++-5

$(NM)/Makefile: $(NM_VER)/Makefile
	mv $(NM_VER) $(NM)
#	cd $@ && patch -p1 < ../snap_patch.diff
	cd $(NM) && patch -p1 < ../snap_features_patch.diff

$(NM_VER)/Makefile: $(AR)
	tar xvfz $<

$(AR):
	wget -O $(AR) $(URL)

.PHONY: clean
clean:
	rm -fr snap snap-$(VER)
	rm -f *.tar.gz snap-aligner
