# Local Branch Sync Plan

## Goal

Preserve a long-lived local-running branch while periodically pulling improvements from `main`.

This branch family exists so the model can run on a local laptop with limited resources. It is not intended to be merged back into `main` unless maintainers explicitly decide to upstream part of the local-running behavior.

## Branch Strategy

Do not create a clean branch from `origin/main` for this purpose. That strategy is only appropriate for upstreaming small isolated changes, such as documentation, into `main`.

For local-running sync work:

1. Start from the current local-running branch, such as `codex-onboarding` or `runlocal`.
2. Create a temporary integration branch.
3. Merge `origin/main` into that integration branch.
4. Resolve conflicts in favor of preserving local laptop execution while adopting main's structural cleanup where possible.
5. Validate the local-running workflow.
6. If successful, keep the integration branch as the refreshed local-running branch or fast-forward/replace the previous local-running branch intentionally.

Suggested commands:

```bash
git fetch origin
git switch codex-onboarding
git switch -c runlocal-main-sync
git merge origin/main
```

If the local-running base should instead be `runlocal`, switch to that branch before creating the integration branch.

## Expected Conflict Areas

A dry merge check from `codex-onboarding` to current `origin/main` predicted conflicts in:

- `fit_objects.Rds`
- `updated_workflow.Rmd`

Main also made broad structural changes, including additions under `core_data/`, `scripts/`, and new helper modules under `R/`. Even if those files do not conflict mechanically, they should be reviewed for overlap with local-running behavior.

## Resolution Principles

- Preserve the branch's ability to run locally with limited cores, memory, and generated-output size.
- Prefer main's cleanup for generated data placement, deleted artifacts, and reorganized shared utilities unless that breaks local execution.
- Do not resurrect large generated artifacts just because they exist on the local-running branch.
- Keep branch-specific documentation clear that this branch is not intended to merge into `main`.
- Avoid mixing scientific/model changes with merge mechanics. If a model assumption must change, make it explicit in a separate commit.

## File-Specific Guidance

### `fit_objects.Rds`

This is a generated/scientific artifact. If `main` deletes or relocates it, prefer accepting main's deletion unless the local workflow cannot regenerate or locate an equivalent input.

If local execution still needs a small laptop-friendly `fit_objects.Rds`, document why it is retained and whether it is a fixture, a cached artifact, or an expected user-generated file.

### `updated_workflow.Rmd`

This is the main manual conflict area.

Preserve local-running affordances such as:

- low-resource settings,
- local batch execution,
- reduced replicate counts or population caps when intentionally used,
- paths that work on a laptop,
- instructions that do not require SLURM.

Adopt from `main` where possible:

- cleaned data locations,
- shared helper functions,
- updated object schemas,
- bug fixes in lineage mapping or karyotype preprocessing,
- documentation comments that clarify the canonical workflow.

### `R/run_local.R`

If `main` deletes this file but does not replace local batch behavior, keep or recreate a local runner. If `main` introduces a better script under `scripts/`, consider wrapping or documenting that instead of maintaining duplicate logic.

### `R/estimate_pmis.R` and `ksim2.cpp`

Review these carefully because they control scientific behavior and performance. Prefer main's bug fixes and simulator improvements unless they make laptop execution infeasible.

If local changes reduce population size, replicate count, grid size, or recording frequency, label them as resource-tuning choices rather than scientific defaults.

### Documentation

After conflict resolution, update:

- `AGENTS.md`
- `docs/architecture.md`
- this sync plan, if the branch strategy changes

The docs should say which branch is meant for local execution and which workflow entry point should be used.

## Validation Checklist

Before declaring the sync successful:

1. Confirm the worktree contains only intentional changes.
2. Confirm local-run commands are documented and still point at existing files.
3. Run a minimal smoke test if feasible, such as sourcing `R/estimate_pmis.R` and running a tiny simulator call.
4. Run the smallest practical local workflow slice, not the full production simulation.
5. Confirm generated outputs either remain untracked or are intentionally documented as branch artifacts.
6. Check that no credentials, `.Rhistory`, `.DS_Store`, logs, or large accidental outputs are staged.

## Open Questions

- Should the refreshed local-running branch be named `runlocal-main-sync`, `runlocal`, or something else after validation?
- Should local-running resource settings be command-line parameters instead of branch-specific code differences?
- Should local fixture data live in git, or should the branch require users to generate it locally from `core_data/`?
- Which outputs are acceptable to track on the local-running branch, if any?

