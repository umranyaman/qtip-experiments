#!/bin/sh

URL=http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/rmsk.txt.gz

for CHR in 1 6 19 22 ; do
    curl ${URL} | gzip -dc | \
	awk "\$6 == \"chr${CHR}\" && (\$12 == \"Simple_repeat\" || \$12 == \"Satellite\" || \$12 == \"Low_complexity\") {print ${CHR},\$7,\$8}" > \
	    rmsk_${CHR}.bed
done

