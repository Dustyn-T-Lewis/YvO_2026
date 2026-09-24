# YvO 2026 · skeletal muscle proteomics, younger against older adults

Analysis code for the resistance-training proteomics study. Vastus lateralis
biopsies from 32 participants, paired pre and post, quantified by DIA mass
spectrometry and analysed with limma.

    Rscript setup.R      # once after cloning
    Rscript run_all.R    # 22 steps, each in its own R session

The analysed matrix is 2,106 proteins across 62 samples, from 32 participants less
the two consensus outliers `Y_S05_Pre` and `Y_S07_Post`. Training changed 135
proteins at FDR < 0.05 in younger adults and none in older adults. Every
headline count, with the script that produces it, is in
`03_DEP/README.md`.

## Layout

```
00_input/          source data; no script writes here
01_normalization/  filter, outlier consensus, cyclic loess   -> S8 Table
02_imputation/     missForest, with the benchmark that chose it -> S9 Table
03_DEP/            limma, four contrasts, and the supplement layer -> S10 Table
04_Figures/
  F00              pipeline QC, supplementary S1a and S1b
  F01 .. F06       Figures 1 to 6; the directory number is the figure number
  abstract_panels  the graphical abstract
  shared/          palettes, helpers, caches; sourced, never run
```

Every stage and figure directory carries its own `README.md` naming what it
reads, what it writes, and which orderings are load-bearing. Read those before
changing anything: several stages depend on files written by an earlier step in
the same run.

Each directory follows one shape. `a_script/` holds the code, `c_data/`
the tables it writes, `b_reports/` the rendered figures. Git tracks the
figures the manuscript cites and the workbooks; every other render regenerates
from a run and is ignored.

## Reproducing

`setup.R` installs from CRAN, Bioconductor and two GitHub sources, then
`renv::restore()` pins every version from `renv.lock`. `run_all.R` runs each
step as a separate `Rscript` child and logs timings to `.runlogs/`.

A full run takes about 16 minutes. The two slowest steps are the effect-size
bootstrap in `03_DEP/a_script/supp/01` at 160 s and the leave-one-subject-out network
refit behind Figure 6 at 255 s.

Scripts prefixed with `_` are sourced by their stitcher and will fail if run on
their own.

## Citation

Archived at [10.5281/zenodo.19886624](https://doi.org/10.5281/zenodo.19886624).
Code is MIT licensed; see `LICENSE`.
