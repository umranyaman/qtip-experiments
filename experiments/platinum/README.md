### Tabulating ROCs

Once we've generated all the ROCs, we want to put them into forms that are easy to summarize into tables and eventually plots.  First let's examine what tables we're targeting.

#### F scores

We want a table that shows, for a range of Fβ's, where β varies the relative importance of precision in the precision/recall calculation.

| β    | Best F<sub>β</sub> original (MAPQ cutoff, QUAL cutoff) | Best F<sub>β</sub> Qtip (MAPQ cutoff, QUAL cutoff) | Relative change |
|------|------------|------------|------------|
| 0.5  |       0.01 () |       0.16 |     + 2.83 |
| 0.75 |       0.26 () |       5.38 |    + 94.37 |
| 1.00 |       0.26 () |       5.38 |    + 94.37 |
| 1.50 |       0.26 () |       5.38 |    + 94.37 |
| 2.00 |       0.26 () |       5.38 |    + 94.37 |
