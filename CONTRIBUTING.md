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

1. Bump `VERSION` on `develop` (e.g. `v0.2.0`). The `update-version.yaml`
   workflow propagates the new version into `README.md`, `pyproject.toml`,
   `src/shiba2corr/__init__.py`, `docker/Dockerfile`, and the issue
   templates, and commits the result back to `develop`.
2. Open a PR from `develop` into `main` with a release-note body.
3. Merge the PR. The `release.yaml` workflow tags the release, creates a
   GitHub Release, and builds & pushes the Docker image. The `publish.yaml`
   workflow uploads the wheel to PyPI.

## Code style

- Tabs for indentation (consistent with the rest of the codebase).
- Keep public-facing APIs in `src/shiba2corr/__init__.py` stable across
  patch releases.
- Add a `pytest` test for any new behavior in `tests/`.

## Reporting issues

Please use the [issue templates](.github/ISSUE_TEMPLATE/) and run
`shiba2corr --verbose` when reporting bugs so we get a full log.
