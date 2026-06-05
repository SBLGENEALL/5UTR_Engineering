# Versioning and Branch Policy

## Official Rule

This repository should be managed with a simple research-pipeline release model.

```text
main        = latest validated working version
v1.0, v1.1  = version tags/releases, not long-lived branches
work/*      = temporary development branches
codex/*     = temporary AI-assisted branches
feature/*   = temporary development branches
```

## Branch Policy

Keep:

- `main`
- temporary active development branches only while they are being worked on

Delete after merge/validation:

- `pr-*`
- `codex/*`
- `feature/*`
- `work/*`
- old experiment branches that are already represented in `main` or a version tag

## Version Policy

Use tags/releases for validated snapshots:

```text
v1.0 = first validated release
v1.1 = validated update
v1.2 = next validated update
```

Do not keep `v1.1`, `v1.2`, etc. as long-lived branches unless there is a specific maintenance reason.

## Pull Request Policy

PRs are temporary review records.

After a PR is merged or rejected:

- keep the PR record closed/merged for history
- delete the source branch if it is no longer active
- record validated changes in `CHANGELOG.md`
- update `MASTER.md` only when the project consensus changes

## Practical Workflow

```text
1. Create temporary branch
2. Develop/test
3. Merge into main after validation
4. Update MASTER.md / CHANGELOG.md / NEXT_ACTIONS.md
5. Create tag/release if this is a validated milestone
6. Delete temporary branch
```
