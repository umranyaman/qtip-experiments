$1 !~ /threshold/ && $2 > 0 && $2 > $3 {
    # SNVs
    tp=$2-$3;
    fp=$3;
    fn=$4;
    snv_prec=tp/(tp+fp);
    snv_recall=tp/(tp+fn);
    snv_f025 = (1 + 0.25 * 0.25) * (snv_prec * snv_recall) / ((0.25 * 0.25) * snv_prec + snv_recall);
    snv_f033 = (1 + (1/3.0) * (1/3.0)) * (snv_prec * snv_recall) / (((1/3.0) * (1/3.0)) * snv_prec + snv_recall);
    snv_f050 = (1 + 0.5 * 0.5) * (snv_prec * snv_recall) / ((0.5 * 0.5) * snv_prec + snv_recall);
    snv_f100 = (1 + 1.0 * 1.0) * (snv_prec * snv_recall) / ((1.0 * 1.0) * snv_prec + snv_recall);
    snv_f200 = (1 + 2.0 * 2.0) * (snv_prec * snv_recall) / ((2.0 * 2.0) * snv_prec + snv_recall);
    snv_f300 = (1 + 3.0 * 3.0) * (snv_prec * snv_recall) / ((3.0 * 3.0) * snv_prec + snv_recall);
    snv_f400 = (1 + 4.0 * 4.0) * (snv_prec * snv_recall) / ((4.0 * 4.0) * snv_prec + snv_recall);
    if(snv_f025 > snv_f025max) {snv_f025max = snv_f025; snv_f025tp = tp; snv_f025fp = fp; snv_f025fn = fn; };
    if(snv_f033 > snv_f033max) {snv_f033max = snv_f033; snv_f033tp = tp; snv_f033fp = fp; snv_f033fn = fn; };
    if(snv_f050 > snv_f050max) {snv_f050max = snv_f050; snv_f050tp = tp; snv_f050fp = fp; snv_f050fn = fn; };
    if(snv_f100 > snv_f100max) {snv_f100max = snv_f100; snv_f100tp = tp; snv_f100fp = fp; snv_f100fn = fn; };
    if(snv_f200 > snv_f200max) {snv_f200max = snv_f200; snv_f200tp = tp; snv_f200fp = fp; snv_f200fn = fn; };
    if(snv_f300 > snv_f300max) {snv_f300max = snv_f300; snv_f300tp = tp; snv_f300fp = fp; snv_f300fn = fn; };
    if(snv_f400 > snv_f400max) {snv_f400max = snv_f400; snv_f400tp = tp; snv_f400fp = fp; snv_f400fn = fn; };

    # Indels
    tp=$5-$6;
    fp=$6;
    fn=$7;
    ind_prec=tp/(tp+fp);
    ind_recall=tp/(tp+fn);
    ind_f025 = (1 + 0.25 * 0.25) * (ind_prec * ind_recall) / ((0.25 * 0.25) * ind_prec + ind_recall);
    ind_f033 = (1 + (1/3.0) * (1/3.0)) * (ind_prec * ind_recall) / (((1/3.0) * (1/3.0)) * ind_prec + ind_recall);
    ind_f050 = (1 + 0.5 * 0.5) * (ind_prec * ind_recall) / ((0.5 * 0.5) * ind_prec + ind_recall);
    ind_f100 = (1 + 1.0 * 1.0) * (ind_prec * ind_recall) / ((1.0 * 1.0) * ind_prec + ind_recall);
    ind_f200 = (1 + 2.0 * 2.0) * (ind_prec * ind_recall) / ((2.0 * 2.0) * ind_prec + ind_recall);
    ind_f300 = (1 + 3.0 * 3.0) * (ind_prec * ind_recall) / ((3.0 * 3.0) * ind_prec + ind_recall);
    ind_f400 = (1 + 4.0 * 4.0) * (ind_prec * ind_recall) / ((4.0 * 4.0) * ind_prec + ind_recall);
    if(ind_f025 > ind_f025max) {ind_f025max = ind_f025; ind_f025tp = tp; ind_f025fp = fp; ind_f025fn = fn; };
    if(ind_f033 > ind_f033max) {ind_f033max = ind_f033; ind_f033tp = tp; ind_f033fp = fp; ind_f033fn = fn; };
    if(ind_f050 > ind_f050max) {ind_f050max = ind_f050; ind_f050tp = tp; ind_f050fp = fp; ind_f050fn = fn; };
    if(ind_f100 > ind_f100max) {ind_f100max = ind_f100; ind_f100tp = tp; ind_f100fp = fp; ind_f100fn = fn; };
    if(ind_f200 > ind_f200max) {ind_f200max = ind_f200; ind_f200tp = tp; ind_f200fp = fp; ind_f200fn = fn; };
    if(ind_f300 > ind_f300max) {ind_f300max = ind_f300; ind_f300tp = tp; ind_f300fp = fp; ind_f300fn = fn; };
    if(ind_f400 > ind_f400max) {ind_f400max = ind_f400; ind_f400tp = tp; ind_f400fp = fp; ind_f400fn = fn; };
}
END {
    print "0.25",snv_f025max,snv_f025tp,snv_f025fp,snv_f025fn,ind_f025max,ind_f025tp,ind_f025fp,ind_f025fn
    print "0.33",snv_f033max,snv_f033tp,snv_f033fp,snv_f033fn,ind_f033max,ind_f033tp,ind_f033fp,ind_f033fn
    print "0.50",snv_f050max,snv_f050tp,snv_f050fp,snv_f050fn,ind_f050max,ind_f050tp,ind_f050fp,ind_f050fn
    print "1.00",snv_f100max,snv_f100tp,snv_f100fp,snv_f100fn,ind_f100max,ind_f100tp,ind_f100fp,ind_f100fn
    print "2.00",snv_f200max,snv_f200tp,snv_f200fp,snv_f200fn,ind_f200max,ind_f200tp,ind_f200fp,ind_f200fn
    print "3.00",snv_f300max,snv_f300tp,snv_f300fp,snv_f300fn,ind_f300max,ind_f300tp,ind_f300fp,ind_f300fn
    print "4.00",snv_f400max,snv_f400tp,snv_f400fp,snv_f400fn,ind_f400max,ind_f400tp,ind_f400fp,ind_f400fn
}
