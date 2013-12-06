#
# wgsim simulator
#
# BTL: I patched it so that it additionally reports whether the mates
# are flipped (i.e. whether mate 2 is to the left of mate 1 w/r/t the
# Watson strand)
#

GIT=$(shell which git)

wgsim: wgsim-git/wgsim
	ln -s -f wgsim-git/wgsim wgsim

wgsim-git/wgsim: wgsim-git
	cd wgsim-git && gcc -O2 -Wall -o wgsim wgsim.c -lz -lm

wgsim-git: $(GIT)
	$(GIT) clone git@github.com:lh3/wgsim.git
	mv wgsim $@
	cd $@ && patch -p1 < ../wgsim.patch
