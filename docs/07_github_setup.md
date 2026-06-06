# GitHub Repository Setup

## Option 1: Create On GitHub First

1. Open GitHub and create a new empty repository named `microc2tracks`.
2. Do not initialize it with a README if you are uploading this repository as-is.
3. On the server or local machine:

```bash
cd microc2tracks
git init
git add .
git commit -m "Initial microc2tracks workflow"
git branch -M main
git remote add origin https://github.com/MichalGd/microc2tracks.git
git push -u origin main
```

## Option 2: Using GitHub CLI

```bash
cd microc2tracks
git init
git add .
git commit -m "Initial microc2tracks workflow"
gh repo create MichalGd/microc2tracks --public --source=. --remote=origin --push
```

Use `--private` instead of `--public` if the repository should not be public yet.

## Upload From The Packaged Archive

```bash
tar -xzf microc2tracks_20260606.tar.gz
cd microc2tracks
git init
git add .
git commit -m "Initial microc2tracks workflow"
git branch -M main
git remote add origin https://github.com/MichalGd/microc2tracks.git
git push -u origin main
```

## Recommended First Checks After Clone

```bash
bash tests/check_bash_syntax.sh
python -m py_compile scripts/validate_samplesheet.py scripts/compare_matrices.py scripts/plot_ps_curves.py
python scripts/validate_samplesheet.py config/samplesheet_example.csv
```

## Linux Line Endings

The repository includes `.gitattributes` so shell scripts are stored with LF line endings. This matters because scripts edited on Windows can otherwise fail on Linux with `bad interpreter` or `$'\r': command not found`.

