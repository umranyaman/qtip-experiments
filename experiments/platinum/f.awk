$1 !~ /threshold/ && $2 > 0 && $2 > $3 {
    tp=$2-$3;
    fp=$3;
    fn=$4;
    prec=tp/(tp+fp);
    recall=tp/(tp+fn);
    f050 = (1 + 0.5 * 0.5) * (prec * recall) / ((0.5 * 0.5) * prec + recall);
    f100 = (1 + 1.0 * 1.0) * (prec * recall) / ((1.0 * 1.0) * prec + recall);
    f200 = (1 + 2.0 * 2.0) * (prec * recall) / ((2.0 * 2.0) * prec + recall);
    if(f050 > f050max) {f050max = f050; f050tp = tp; f050fp = fp; f050fn = fn; };
    if(f100 > f100max) {f100max = f100; f100tp = tp; f100fp = fp; f100fn = fn; };
    if(f200 > f200max) {f200max = f200; f200tp = tp; f200fp = fp; f200fn = fn; };
}
END {
    print "0.5",f050max,f050tp,f050fp,f050fn
    print "1.0",f100max,f100tp,f100fp,f100fn
    print "2.0",f200max,f200tp,f200fp,f200fn
}
