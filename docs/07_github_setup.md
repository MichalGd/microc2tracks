# Manual GitHub Update Procedure

These are operator instructions only. The candidate preparation process does not
push, create a branch, open a pull request, or modify the remote repository.

1. Verify the archive checksum against `SHA256SUMS.txt` and inspect
   `UPLOAD_MANIFEST.txt`.
2. Clone the intended repository into a new directory. Do not reuse the local
   baseline whose configured origin contains the `micorc2tracks` typo.
3. Confirm the remote and starting commit:

```bash
git clone https://github.com/MichalGd/microc2tracks.git microc2tracks-review
cd microc2tracks-review
git remote -v
git status --short
git log -1 --oneline
```

4. Create a review branch locally if desired, then copy only paths listed in
   `UPLOAD_MANIFEST.txt` from the extracted candidate. Do not copy `.git`,
   caches, the external archive/checksum file, or the external comparison patch.
5. Review and test before committing:

```bash
git status --short
git diff --check
git diff --stat
bash tests/check_bash_syntax.sh
bash tests/test_reference_sentinel.sh
python -m py_compile scripts/*.py tests/test_multireference.py
python -m unittest discover -s tests -p 'test_*.py' -v
python scripts/reference_registry.py validate --registry config/references.tsv
python scripts/validate_samplesheet.py \
  --reference-registry config/references.tsv \
  config/samplesheet_example.csv
```

6. Review `config/references.tsv` carefully: its `/shared/references/...` paths
   are public examples, not proof of an installed server layout.
7. Commit with a focused message, inspect the commit, and only then push under
   the repository owner's normal review policy:

```bash
git add -- $(cat UPLOAD_MANIFEST.txt)
git commit -m "Add per-sample mm39 and hg38 reference support"
git show --stat --oneline HEAD
# Push only after explicit authorization and normal project review.
```

The repository includes `.gitattributes` so shell scripts remain LF-terminated
when reviewed on Windows.
