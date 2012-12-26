Bowtie 2 was augmented in a few specific ways to enable these analyses:
- Added the --mapq-extra option, which causes Bowtie 2 to output some
  additional SAM extra fields:
  + ZP:i: Score of best concordant paired-end alignment
  + Zp:i: Second-best concordant paired-end alignment score
  + ZU:i: Score of best unpaired alignment
  + Zu:i: Score of second-best unpaired alignment
  + XP:B:I: String describing seed hits
  + Xs:i: Best invalid alignment score of this mate
  + Ys:i: Best invalid alignment score of opposite mate
- 