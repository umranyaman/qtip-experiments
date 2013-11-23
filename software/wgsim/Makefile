#
# wgsim simulator
#

GIT=$(shell which git)

wgsim: wgsim-git/wgsim
	ln -s -f wgsim-git/wgsim wgsim

wgsim-git/wgsim: wgsim-git
	cd wgsim-git && gcc -O2 -Wall -o wgsim wgsim.c -lz -lm

wgsim-git: $(GIT)
	$(GIT) clone git@github.com:lh3/wgsim.git
	mv wgsim $@
