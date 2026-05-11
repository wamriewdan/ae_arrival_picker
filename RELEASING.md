# Releasing `ae-picker`

This project is configured to publish to PyPI with GitHub Actions Trusted
Publishing.

The repository also has a separate CI workflow that runs on every push and
pull request to validate installation, CLI wiring, and package metadata
without publishing anything.

Relevant workflows:

- `.github/workflows/ci.yml`
- `.github/workflows/publish-testpypi.yml`
- `.github/workflows/publish-pypi.yml`

## One-time setup

1. Confirm the distribution name you want on PyPI is `ae-picker`.
2. Push the workflow files `.github/workflows/publish-pypi.yml` and
   `.github/workflows/publish-testpypi.yml` to GitHub.
3. In GitHub, create these environments for this repository:
   - `testpypi`
   - `pypi`
4. In TestPyPI, configure Trusted Publishing for this exact workflow:
   - TestPyPI project name: `ae-picker`
   - Repository owner: `wamriewdan`
   - Repository name: `ae_arrival_picker`
   - Workflow file name: `publish-testpypi.yml`
   - Environment name: `testpypi`
5. In PyPI, configure Trusted Publishing for this exact workflow:
   - PyPI project name: `ae-picker`
   - Repository owner: `wamriewdan`
   - Repository name: `ae_arrival_picker`
   - Workflow file name: `publish-pypi.yml`
   - Environment name: `pypi`
6. If the project does not exist on TestPyPI or PyPI yet, create a pending
   publisher from the account publishing settings page on that index instead of
   the project settings page.

## TestPyPI dry run

1. Make sure `pyproject.toml` contains the version you want to test.
2. Commit and push the version change to `main`.
3. In GitHub, run the `Publish to TestPyPI` workflow manually.
4. Wait for the workflow to finish successfully.
5. Verify the TestPyPI project page.
6. Verify installation in a clean environment:

```bash
python -m venv .venv-testpypi
. .venv-testpypi/Scripts/Activate.ps1
python -m pip install -U pip
python -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple ae-picker
python -c "import ae_picker; print(ae_picker.__version__)"
ae-picker --help
```

`--extra-index-url https://pypi.org/simple` is useful because your package
will come from TestPyPI but most dependencies are usually resolved from PyPI.

## Release checklist

1. Make sure `pyproject.toml` has the version you want to release.
2. Review `README.md` so the install and usage examples match the release.
3. Optionally run the `Publish to TestPyPI` workflow first.
4. Commit and push the release changes to `main`.
5. Create a Git tag that matches the version, for example `v0.2.0`.
6. Create a GitHub Release from that tag and publish it.
7. Wait for the `Publish to PyPI` workflow to finish successfully.
8. Verify the package page on PyPI.
9. Verify installation in a clean environment:

```bash
python -m venv .venv-test
. .venv-test/Scripts/Activate.ps1
python -m pip install -U pip
python -m pip install ae-picker
python -c "import ae_picker; print(ae_picker.__version__)"
ae-picker --help
```

## First release note

For the first PyPI release, Trusted Publishing does not reserve the package
name ahead of time. A pending publisher only becomes active when the first
publish succeeds.

## If a release fails

1. Open the failed GitHub Actions run and inspect the failing step.
2. If PyPI says the publisher is invalid, double-check:
   - repository owner
   - repository name
   - workflow filename
   - environment name
3. If the upload failed after files were accepted, do not rebuild the same
   version with different contents. Bump the version and release again.
