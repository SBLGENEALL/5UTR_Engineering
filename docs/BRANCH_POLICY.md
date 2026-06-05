# Branch Policy

Last updated: 2026-06-04

## Current branch roles

```text
main
= stable release branch
= do not use for active PR experiments

improved-v1.2
= current validated working baseline
= PR2 validation passed
= recommended base for PR3

improved-v1.2-pr3-selection-policy
= next active PR3 development branch
= selection policy / uAUG audit / dry-run work
```

## Cleanup rule

PR branches are temporary workspaces.

After a PR branch is validated and promoted to a versioned baseline branch:

1. Record the validation result in `docs/MASTER.md`.
2. Record the change and validation summary in `docs/CHANGELOG.md`.
3. Create or update the version baseline branch, for example `improved-v1.2`.
4. Delete the old PR branch.

## Branches that can be deleted after v1.2 promotion

The following PR/work branches are considered temporary and can be removed once `improved-v1.2` is confirmed present:

```text
improved-v1.1-pr2-tss-expression
improved-v1.1-pr1.1
```

Keep:

```text
main
improved-v1.1
improved-v1.2
improved-v1.2-pr3-selection-policy
```

Optional: `improved-v1.1` may be kept as the historical PR1/PR1.1 validated baseline. If branch count becomes too high later, it can be archived after `improved-v1.2` and later versions are stable.

## Manual cleanup commands

From a local clone with GitHub push permission:

```bash
git fetch --all --prune

git push origin --delete improved-v1.1-pr2-tss-expression
git push origin --delete improved-v1.1-pr1.1

git branch -r
```

Do not delete `main` or the current validated version branch.
