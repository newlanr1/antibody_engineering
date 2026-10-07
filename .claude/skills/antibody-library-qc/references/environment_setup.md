# Environment setup for Analyses 5 & 6

Sapiens (humanness) and IgFold (structure) both need Python versions incompatible
with this skill's main venv (3.9.6) *and* with each other — Sapiens needs <=3.8,
IgFold needs >=3.11. Each lives in its own dedicated environment that the main venv
calls into via `subprocess`, since you can't import across Python versions directly.
**Both bridge scripts (`humanness_scoring.py`, `structure_prediction.py`) skip
gracefully with a clear message if their environment isn't set up** — this is a
one-time, manual, per-machine prerequisite, not something every run can assume.

**IgFold (Python 3.11) — low risk, plain Homebrew + pip:**
```bash
brew install python@3.11
/opt/homebrew/bin/python3.11 -m venv venv-igfold311
source venv-igfold311/bin/activate && pip install igfold pandas
```

**Sapiens (Python 3.8) — not straightforward.** Building CPython 3.8 from source
(`pyenv install 3.8.20`) segfaults on modern macOS/Apple Silicon — confirmed on this
machine; the compiled interpreter crashes just being invoked, not merely during the
pip-bootstrap step. `conda`/`miniforge` was tried as a fallback and hit its own
problem: conda's bundled CA store doesn't trust this network's TLS-intercepting
corporate proxy on `anaconda.org`/`prefix.dev` (while PyPI/GitHub/Hugging
Face/Zenodo aren't intercepted, so plain `pip` against those domains is unaffected).
What actually worked: a **prebuilt** (not locally-compiled) CPython distribution from
[`astral-sh/python-build-standalone`](https://github.com/astral-sh/python-build-standalone),
which sidesteps both problems since nothing is compiled locally and nothing touches
conda's network stack:
```bash
curl -L -o /tmp/cpython38.tar.gz \
  https://github.com/astral-sh/python-build-standalone/releases/download/20231002/cpython-3.8.18+20231002-aarch64-apple-darwin-install_only.tar.gz
tar -xzf /tmp/cpython38.tar.gz -C /tmp
mkdir venv-sapiens38 && cp -R /tmp/python/* venv-sapiens38/
venv-sapiens38/bin/python3.8 -m pip install sapiens pandas
```
(That specific release tag is one confirmed to still carry a 3.8 build — astral-sh
drops old versions over time, so check their releases page for a current tag with an
`aarch64-apple-darwin` 3.8 asset if this one 404s. For Intel Macs, swap in the
`x86_64-apple-darwin` asset name instead.) Despite its README saying 3.9+ breaks it
"due to a fairseq bug," the pip-installed `sapiens` package itself has no `fairseq`
dependency and worked fine once a real 3.8 interpreter was available — the README's
caveat seems to be about an unrelated dev/training path, not the inference API used
here.

Both venv paths are auto-detected by `humanness_scoring.py` / `structure_prediction.py`
relative to the repo root (`venv-sapiens38/bin/python3.8`, `venv-igfold311/bin/python`)
— override with `--sapiens-python`/`--igfold-python` if yours live elsewhere. Both
directories are already in `.gitignore`.

## If you're setting this up on a different machine

The specific failures above (CPython 3.8 segfaulting on build, conda's CA store not
trusting a corporate TLS proxy) are environment-specific — a different machine might
hit neither and `pyenv install 3.8.20` / `conda create -n sapiens38 python=3.8` might
just work. Try the straightforward path first; fall back to the
python-build-standalone approach only if you hit a similar wall.
