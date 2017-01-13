#!gmake

NM=varmatch
USER=medvedevgroup
EXE=varmatch

.PHONY: all
all: $(NM)

$(EXE): $(NM)-git/$(EXE)
	ln -s -f $< $@

$(NM)-git/$(EXE): $(NM)-git
	$(MAKE) -C $<

$(NM)-git:
	git clone --recursive https://github.com/$(USER)/$(NM).git $@

.PHONY: clean
clean:
	rm -rf $(NM)-git $(NM)
