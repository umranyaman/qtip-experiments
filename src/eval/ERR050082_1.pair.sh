#!/bin/sh

cat >.ERR050082_1.pair.sh <<EOF
#!/bin/sh
make ERR050082_1.pair.csv
EOF

sbatch .ERR050082_1.pair.sh -N 1-1 -n 16 --time=48:00:00 --mem=32g -p shared
