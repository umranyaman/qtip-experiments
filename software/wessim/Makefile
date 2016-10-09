#!gmake

NM=Wessim
AR=$(NM).tar.gz
URL=https://github.com/sak042/$(NM)/archive/$(AR)

$(NM): $(AR)
	tar jxvf $<
	mv $(NM)-$(NM) $@

$(AR):
	wget $(URL)

.PHONY: clean
clean:
	rm -fr $(NM)
	rm -f $(AR)
