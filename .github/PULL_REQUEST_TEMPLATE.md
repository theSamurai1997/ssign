# Pull request

## Summary

<!-- One or two sentences: what does this PR change and why? -->

## Linked issue

<!-- e.g. Closes #123, Related to #456. Open an issue first for non-trivial changes. -->

## Type of change

- [ ] Bug fix (no behaviour change for users without the bug)
- [ ] Feature (adds new functionality)
- [ ] Refactor (no functional change)
- [ ] Documentation
- [ ] CI / build / packaging

## How was this tested?

<!--
Required for code changes. Note specific test files added or modified, and
any manual verification steps. For pipeline changes, mention which input
genomes you ran against and what changed in the output.
-->

## Checklist

- [ ] Code follows the conventions in `CONTRIBUTING.md`
- [ ] New unit tests added (or existing coverage extended)
- [ ] `pytest tests/unit/` passes locally
- [ ] `ruff check src/ tests/` passes locally
- [ ] `mypy` passes locally
- [ ] `CHANGELOG.md` updated under `[Unreleased]` for user-facing changes
- [ ] Documentation updated where relevant
- [ ] If touching `docs/`, follows the Diátaxis section (tutorial / how-to / reference / explanation)
