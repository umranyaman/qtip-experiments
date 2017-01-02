#!gmake

NM=vcflib
USER=ekg
VER=1.0.0-rc1
EXE=vcfroc

.PHONY: all
all: $(NM)

$(EXE): $(NM)-git/bin/$(EXE)
	ln -s -f $< $@

$(NM)-git/bin/$(EXE): $(NM)-git
	$(MAKE) -C $<

$(NM)-git:
	git clone -b v$(VER) --recursive https://github.com/$(USER)/$(NM).git $@

.PHONY: clean
clean:
	rm -rf $(NM)-git $(NM)
