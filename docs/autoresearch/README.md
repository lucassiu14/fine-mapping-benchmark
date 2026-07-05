# `docs/autoresearch/` — iteration log for the auto-research loop

This directory hosts the running record of the auto-research loop
(fmbenchmark_autoresearch_updated). One markdown file per iteration:

- `iteration-000-phase0.md` — Phase 0 setup + method additions (this loop's
  initial pass).
- `iteration-001-phase1-grid.md` — Phase 1 HPC grid + worker plumbing.
- `iteration-NNN-<slug>.md` — each subsequent Phase 2 iteration.

Novel methods get their own descriptor file:

- `method-<method_name>.md`

Conventions per §2.6 of the plan:

- Iteration number in the filename and body header.
- Issues encountered / steps taken.
- Methods tested (with type tags: `baseline`, `hyperparameter_variant`,
  `novel`) and non-default arg info.
- Key insights about the methods AND the datasets.
- Names of new methods added to the feasible set for the next iteration.
- For novel methods designed but not yet built: proposal in the log,
  followed by a `method-<name>.md` file once implemented.

Memory update: before every break, update the persistent memory (see the
`memory/` notes in the parent project) with where we are, results
we're waiting on, and specific next steps — so a fresh session can
resume purely from the on-disk state.
