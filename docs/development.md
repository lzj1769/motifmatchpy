# Contributing and publishing docs

## Develop the package

```sh
git clone https://github.com/lzj1769/motifmatchpy.git
cd motifmatchpy
uv sync
uv run pytest
uv run pytest -m reference
uvx ruff@0.16.6 check src tests
```

`uv sync` builds the development-only C++ reference, so a C++17 compiler is
required. Runtime users do not need that reference. Tests include independent
brute-force checks and comparisons against the patched reference. Use four-space
indentation, add regression tests for algorithm changes, and describe the change
and validation in pull requests. Record reference patches in
`reference/core/VENDORED.md`.

## Build or preview documentation

From the repository root, create a separate environment:

```sh
python -m venv .venv-docs
.venv-docs/bin/python -m pip install -r docs/requirements.txt .
.venv-docs/bin/python scripts/check_docs_examples.py
.venv-docs/bin/python -m mkdocs build --strict
.venv-docs/bin/python -m mkdocs serve
```

On Windows, use `.venv-docs\Scripts\python` instead. The preview is available at
`http://127.0.0.1:8000/`. Built HTML goes to `site/`, which is ignored by Git.
API reference pages are generated from `src/` by mkdocstrings. After changing
runtime code, reinstall the package in the docs environment before running
example checks so they use the same implementation as the reference pages.

Edit Markdown under `docs/` and update `mkdocs.yml` when adding navigation pages.
The example checker runs the CLI commands and Python snippets in the guide in
a temporary directory, checks inline assertions, and compares published stdout
examples. Installation commands are intentionally not executed by that checker.

## GitHub Pages workflow

The [Documentation workflow](https://github.com/lzj1769/motifmatchpy/actions/workflows/docs.yml)
builds on documentation or package changes in pull requests and on `main`.
It also supports manual dispatch. Builds run the example checker and
`mkdocs build --strict`; warnings fail the build.

Only successful builds on `main` deploy. The deploy job:

1. Downloads the already-built site artifact.
2. Commits the HTML to the `gh-pages` branch with `ghp-import`, preserving history.
3. Publishes the same site with GitHub's Pages deployment actions.

Pull requests never receive the deploy job's write permissions. No personal
access token is stored in the workflow; deployment uses the repository's
`GITHUB_TOKEN` and Pages OIDC support.

### One-time repository setup

An administrator must select **Settings → Pages → Build and deployment →
Source: GitHub Actions**. Allow the `github-pages` environment to deploy from
`main`, and permit the workflow to write the `gh-pages` branch if branch rules
restrict it. Then push to `main` or run the Documentation workflow manually.

The intended URL is <https://lzj1769.github.io/motifmatchpy/>. The workflow also
maintains `gh-pages` as a copy of the generated site, but publication uses
Actions: pushes made with `GITHUB_TOKEN` do not themselves trigger a branch-based
Pages build. See GitHub's
[publishing-source documentation](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site)
and [custom-workflow guide](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages).

If deployment fails, inspect the deploy job, Pages source setting, environment
rules, and branch restrictions. Do not edit HTML on `gh-pages` manually: the
next successful deployment regenerates it from `docs/`.
