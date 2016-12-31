#!/usr/bin/env bash

echo "Making real-data SAM links..."
for pa in unp pair ; do
    for aln in bt2 bwa snap ; do
        for dat in 2 3 ; do
            FN="ERR05008${dat}_1.${aln}.${pa}.sam"
            if [ ! -f "${FN}" ] ; then
                echo "  * Making link to ../real_data/${FN}..."
                if [ ! -f "../real_data/${FN}" ] ; then
                    echo "  ERROR: ../real_data/${FN} does not exist"
                    exit 1
                fi
                ln -s -f ../real_data/${FN} ${FN}
            else
                echo "  Link to ../real_data/${FN} already present"
            fi
        done
    done
done

echo "Making simulated data links..."
for rdlen in 100 250 ; do
    for genome in hg mm ; do
        for aln in bt2s snap bwamem ; do
            # Unpaired
            FN="r0_${aln}_${genome}_mason_ill_${genome}_${rdlen}"
            if [ ! -f "${FN}.trial0.sam" ] ; then
                echo "  * Making link to ../simulated_reads/various_genomes/${FN}.out/trial0/final.sam..."
                if [ ! -f "../simulated_reads/various_genomes/${FN}.out/trial0/final.sam" ] ; then
                    echo "  ERROR: ../simulated_reads/various_genomes/${FN}.out/trial0/final.sam does not exist"
                    exit 1
                fi
                ln -s -f ../simulated_reads/various_genomes/${FN}.out/trial0/final.sam ${FN}.trial0.sam
            else
                echo "  Link to ../simulated_reads/various_genomes/${FN}.out/trial0/final.sam already present"
            fi

            # Paired
            FN="r12_${aln}${rdlen}_${genome}_mason_ill_${genome}_${rdlen}"
            if [ ! -f "${FN}.trial0.sam" ] ; then
                echo "  * Making link to ../simulated_reads/various_genomes/${FN}.out/trial0/final.sam..."
                if [ ! -f "../simulated_reads/various_genomes/${FN}.out/trial0/final.sam" ] ; then
                    echo "  ERROR: ../simulated_reads/various_genomes/${FN}.out/trial0/final.sam does not exist"
                    exit 1
                fi
                ln -s -f ../simulated_reads/various_genomes/${FN}.out/trial0/final.sam ${FN}.trial0.sam
            else
                echo "  Link to ../simulated_reads/various_genomes/${FN}.out/trial0/final.sam already present"
            fi
        done
    done
done
