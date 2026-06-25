# MERGE CES Parametrisation

This repository contains the R workflow used to estimate nested CES production
parameters for MERGE-ETL. The workflow is deliberately staged so the final
MERGE table is traceable back to input-data diagnostics, broad Stage 1 searches,
targeted Stage 1 refinements, and the final high-resolution Stage 2 estimation.

## Folder And File Map

- `MERGE CES workflow guide.md` is the short operational guide. Start there when
  deciding what to run next.
- `MERGE CES parameters.R` is the main estimator. It runs Stage 1, Stage 2, or a
  small integration test depending on the settings at the top of the file.
- `MERGE CES parameters analysis.R` reads a saved RDS result and produces review
  tables, the scientific report, and publication figures.
- `stage0/MERGE CES parameters stage0 analytics.R` audits the raw macro input
  before CES estimation. It does not estimate CES parameters.
- `MERGE descriptive stats.xlsx`, `MERGE regions proposal.xlsx`, and
  `test/MERGE regions proposal.xlsx` are supporting workbooks.
- `stage 1 rho distribution.png` and `Figures.pptx` are presentation artifacts.

The generated `.csv`, `.rds`, and `.RData` outputs are intentionally ignored by
Git because full estimation outputs can be large and are reproducible from the
scripts and input data.

## Scientific Structure

The production function is a two-level nested CES specification. The lower nest
combines capital and labour into value added. The upper nest combines value
added and energy into gross output. The estimator follows the CES convention
introduced by Arrow, Chenery, Minhas, and Solow (1961), where substitution
elasticity is represented separately from distribution weights.

In this workflow, `sigma` is the substitution elasticity and `rho` is the
curvature parameter used by `micEconCES`:

```text
sigma = 1 / (1 + rho)
rho   = 1 / sigma - 1
```

The scripts define grids in `sigma` space because `sigma` is easier to interpret
economically. The grids are converted to `rho` immediately before estimation.

## Recommended Run Order

1. Run Stage 0 when the macro input data changes.
2. Run the broad Stage 1.1 search if no current-format Stage 1 history exists.
3. Run targeted Stage 1 rounds by changing `stage1` to `stage1.2`, `stage1.3`,
   etc. The main script infers the previous history file automatically.
4. Run Stage 2 from `MERGE CES parameters.R`.
5. Run `MERGE CES parameters analysis.R` for review tables, figures, and the
   scientific report.

The final MERGE-facing output is the Stage 2
`stage2_merge_iam_parameter_table.csv`, mirrored inside the Stage 2 RDS as
`merge_iam_parameter_table`.

## Verification Expectations

After editing the scripts, at minimum parse all R files with `Rscript` or
`R.exe` before committing. For estimator changes, also run a small test mode and
check that the generated RDS object still contains the main tables documented in
`MERGE CES workflow guide.md`.

## Reference

Arrow, K. J., Chenery, H. B., Minhas, B. S., and Solow, R. M. (1961). Capital
labor substitution and economic efficiency. Review of Economics and Statistics,
43(3), 225-250.
