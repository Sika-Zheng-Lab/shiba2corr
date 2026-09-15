# Contributing to shiba2corr

Thanks for taking the time to contribute!

## Development setup

```bash
git clone https://github.com/Sika-Zheng-Lab/shiba2corr.git
cd shiba2corr
pip install -e ".[dev]"
pytest -q
```

## Branching workflow

- `main` — release branch; protected.
- `develop` — integration branch; pull requests target this branch.
- Feature branches: `feature/<short-description>` from `develop`.
- Bug fixes: `fix/<short-description>` from `develop`.

Open a pull request from your feature/fix branch into `develop`. CI runs
tests automatically.

## Releasing

1. Bump `VERSION` on `develop` (e.g. `v1.1.0`). The `update-version.yaml`
   workflow propagates the new version into `README.md`, `pyproject.toml`,
   `src/shiba2corr/__init__.py`, `docker/Dockerfile`, and the issue
   templates, and commits the result back to `develop`.
2. Update the release date and notes in `CHANGELOG.md`, then open a PR from
   `develop` into `main`.
3. Merge the `develop` to `main` PR. The `release.yaml` workflow reads
   `VERSION`, creates the corresponding annotated tag on the merge commit,
   validates all version fields, runs the test matrix, builds and attaches the
   Python distributions to a GitHub Release using the PR body as its release
   notes, publishes them to PyPI using Trusted Publishing, and then publishes
   the versioned and `latest` Docker images. A PR from any branch other than
   `develop` does not create a release.

Pushing an existing release tag manually also starts the workflow. This is
intended only for recovery when the automatic workflow could not run.

The GitHub `pypi` environment must be registered as a PyPI Trusted Publisher.
Docker publishing requires the `DOCKER_HUB_USERNAME` and
`DOCKER_HUB_ACCESS_TOKEN` repository secrets.

## Code style

- Tabs for indentation (consistent with the rest of the codebase).
- Keep public-facing APIs in `src/shiba2corr/__init__.py` stable across
  patch releases.
- Add a `pytest` test for any new behavior in `tests/`.

## Reporting issues

Please use the [issue templates](.github/ISSUE_TEMPLATE/) and run
`shiba2corr --verbose` when reporting bugs so we get a full log.
