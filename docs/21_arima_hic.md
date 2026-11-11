# Arima two-enzyme Hi-C

The native `arima_2enzyme` protocol is a Hi-C subvariant: keep `assay=hic`.
It generates a combined restriction-fragment map for `^GATC` and `G^ANTC`
and uses the existing alignment, deduplication, contact filtering and matrix
stages. No separate preprocessing command is required.

This support is included in **4.1**. The older
**4.0** release has neither this preset nor correct degenerate-`N` motif
handling. Publishing this code does not deploy it or change an installed
release. Follow the [4.1 deployment guide](22_version_4_1_deployment.md) before
selecting the new preset on a server; keep active runs on their existing release.

## Scope and configuration

This profile describes the **two-recognition-pattern Arima protocol**, not all
Arima kits or library generations. Confirm the study's protocol before using
it. Different restriction mixtures require their own validated definitions;
the vendor name alone does not establish the appropriate profile.

For a homogeneous project, retain the existing reference, resource and output
settings and set:

```bash
HIC_PROTOCOL="arima_2enzyme"
HIC_RESTRICTION_ENZYMES=""
```

An empty global enzyme setting lets the profile supply `DpnII+HinfI`. Use the
new release's bundled `config/hic_protocols.json`; an older custom registry
does not gain this profile automatically. Literal notation such as
`(^GATC, G^ANTC)` is explanatory notation, **not valid configuration syntax**.
Explicit patterns are supported using `+` as the separator, as shown below.

Alternatively, select `arima_2enzyme` in each relevant samplesheet row's
`hic_protocol` field. The [example samplesheet](../config/samplesheet_arima_example.csv)
uses an explicit `DpnII+HinfI` row value to avoid an unrelated inherited global
enzyme setting. Replace its paths and reference identifier with your registered
reference and real FASTQs. Micro-C rows must leave the Hi-C-specific fields
empty.

Selection precedence is unchanged:

1. The row's `hic_protocol` overrides global `HIC_PROTOCOL`.
2. The row's `restriction_enzymes` overrides global `HIC_RESTRICTION_ENZYMES`.
3. If neither enzyme field is set, the chosen profile supplies its default.

For example, changing only the profile while retaining a global
`HIC_RESTRICTION_ENZYMES="DpnII"` is a conflict. Clear that global value or set
`DpnII+HinfI` explicitly. The profile rejects a different enzyme combination
or changed motif/cleavage definitions rather than silently applying a different
digest. This preserves the behavior of existing configuration overrides.

## Enter restriction patterns directly

For an explicitly specified digest, the existing generic Hi-C profile accepts
cut patterns directly in the existing enzyme field:

```bash
HIC_PROTOCOL="restriction_hic"
HIC_RESTRICTION_ENZYMES="^GATC+G^ANTC"
```

The same syntax works in the samplesheet's `restriction_enzymes` field. No new
column or separate preprocessing action is required. Registered enzyme names
and patterns can be mixed, for example `DpnII+G^ANTC`. Novel patterns do not
need a registry edit: for example, `C^YCGRG` is a valid palindromic IUPAC pattern
with an explicit cut position. Use only the pattern and cut supported by your
actual library protocol; these syntax examples do not identify a new study's
digest automatically.

Each direct pattern must contain exactly one ASCII caret (`^`). It marks the
plus-strand cut between bases: `^GATC` has offset 0, `G^ANTC` has offset 1, and
`CATG^` has offset 4. Lowercase pattern letters are normalized to uppercase.
Do not add parentheses, commas, internal whitespace, a Unicode lookalike
caret, or regular-expression syntax. A bare unregistered motif such as
`GANTC` is rejected because its cleavage position is unspecified. Separate
tokens with `+`; surrounding whitespace around each token is allowed.

Matching direct patterns reuse registered definitions. For a given profile
with the bundled registry, `^GATC+G^ANTC` and `DpnII+HinfI` resolve to the same
policy and stage signatures. The resolver prefers the selected profile's
preset enzyme names, then matching registry names in sorted order. Explicit
existing aliases, such as `MboI`, retain their existing names. Novel patterns
receive stable motif/offset identifiers; conflicting identifier definitions
are rejected.

For Chen's two-pattern library, `arima_2enzyme` remains the preferred preset
because it verifies the expected definitions. It also accepts the equivalent
explicit pattern pair. A contradictory pattern or altered registered cut is
rejected. Choosing generic `restriction_hic` with the same patterns generates
the same digest but records a different profile label in provenance.

## Recognition, cleavage and reference handling

| Registry name | Recognition sequence | Cleavage notation | Plus-strand offset |
| --- | --- | --- | ---: |
| DpnII | `GATC` | `^GATC` | 0 |
| HinfI | `GANTC` | `G^ANTC` | 1 |

These names label computational recognition definitions. They do not establish
which commercial enzymes were supplied in a proprietary cocktail. Cleavage
positions are supported by the manufacturer's [DpnII](https://www.neb.com/en-us/products/r0543-dpnii)
and [HinfI](https://www.neb.com/en-us/products/r0155-hinfi) specificity information.

Both direct patterns and registry definitions accept the full DNA IUPAC
alphabet:

| Symbol | Matching reference bases |
| --- | --- |
| A | A |
| C | C |
| G | G |
| T | T |
| R | A or G |
| Y | C or T |
| S | C or G |
| W | A or T |
| K | G or T |
| M | A or C |
| B | C, G or T |
| D | A, G or T |
| H | A, C or T |
| V | A, C or G |
| N | A, C, G or T |

Each symbol expands to **known A/C/G/T reference bases only**. An unresolved
reference `N`, `R` or other ambiguity code does not establish a restriction
site. Registry definitions use an uppercase motif and separate numeric
`cut_offset`; direct patterns instead place the caret inside/beside the motif.

The complete motif must equal its degenerate reverse complement. The
complement pairs include R/Y, K/M, B/V and D/H; S, W and N complement
themselves. Cuts must lie within the motif or on either boundary. Thus this
extension does **not** support every restriction enzyme: nonpalindromic
recognition patterns and Type IIS cuts outside the motif require different
handling and remain unsupported. Arbitrary regular expressions are rejected.
The scanner handles overlapping sites and sites crossing input chunk/line
boundaries using the biological motif length.

The workflow automatically combines both enzymes' cut coordinates into one
sorted, unique set per chromosome and generates the fragment BED from the
exact configured FASTA. It uses 0-based interbase cleavage coordinates and
half-open BED intervals; 1-based pairs coordinates are converted before
fragment assignment. The digest cache records the reference identities,
enzyme definitions, coordinate convention and generator version (2).

For mm39 remapping, use the validated mm39 FASTA, BWA index and matching
chromosome sizes. An mm10 site file cannot be reused with mm39. The generated
BED follows microc2tracks' cleavage-coordinate convention; it is not intended
to reproduce Juicer's recognition-position text file byte for byte.

## Chen GSE212455 evidence and preprocessing

For Chen et al., *Science Advances* (2023),
[doi:10.1126/sciadv.abo7605](https://doi.org/10.1126/sciadv.abo7605), the
[GEO study](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212455)
and sample protocols describe an Arima kit/cocktail with Swift Accel-NGS 2S
Plus library preparation. The exact kit catalogue/version is not established
by that description alone.

More specific evidence is embedded in the authors' deposited processed
matrices for [WT GSM6533053](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6533nnn/GSM6533053/suppl/GSM6533053_WT_Martix_HiC_inter_30.hic)
and [Tet2/Tet3 DKO GSM6533054](https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6533nnn/GSM6533054/suppl/GSM6533054_TetDKO_Matrix_HiC_inter_30.hic).
Their inspected header statistics record Juicer version 1.6, `-s Arima`,
an mm10 reference and a restriction-site file named `mm10_Arima.txt`.
The version-matched official
[Juicer 1.6 site generator](https://raw.githubusercontent.com/aidenlab/juicer/1.6/misc/generate_site_positions.py)
defines `Arima` as `GATC` plus `GANTC`. Together these observations support
the two-pattern profile for this study with high confidence. The authors'
actual historical site file was not retrieved, so its contents and exact kit
identity remain unverified.

Selecting this profile adds **no BL-Hi-C bridge-linker processing** and no
guessed fixed five-base/UMI trimming. It uses the project's configured ordinary
fastp processing and existing BWA/pairtools Hi-C alignment/parsing. Any
additional library-specific trimming requires its own evidence and explicit
settings; it is not inferred from the Arima brand name.

## Filtering and comparison with other datasets

The profile selects the digest; it does not override project filtering
parameters. Existing defaults remain MAPQ 30, Hi-C minimum cis distance 0,
same-fragment removal enabled, and adjacent inward-fragment removal disabled.
These are microc2tracks defaults, not a claim of exact reproduction of the
authors' Juicer artifact filtering.

Use consistent reference assembly, chromosome selection, mapping-quality
threshold, duplicate policy, matrix resolutions and normalization where
scientifically appropriate. Retain the protocol-specific digest. Matching
these settings improves comparability but cannot remove differences caused by
restriction chemistry, sequencing depth, cell composition or library design.
Hi-C and Micro-C should not be forced to share an inappropriate fragment or
short-range-contact filter merely to make settings identical.

The resulting `.hic` and `.mcool` outputs support common downstream analysis.
They do not by themselves reproduce a paper's original Juicer results,
HiCCUPS/Arrowhead calls or statistical comparisons. A new mm39/pairtools run
is not expected to reproduce the deposited mm10/Juicer contact counts exactly.

## Validation and existing runs

Run the focused native tests from this checkout:

```bash
python3 -B -m unittest discover -s tests -p test_arima.py
```

These tests use small fixtures; they are not a full production mapping run.
Existing literal-motif transforms are preserved. The profile and enzyme
definitions participate in accepted-contact/refilter signatures, so changing
the digest changes those stages and downstream products without intrinsically
changing ordinary alignment or deduplication settings. A cache generator
version change can regenerate a fragment BED when a digest is needed.

The earlier private Chen `4.0+chen-arimaN1` runtime and project-specific registry
are **not migrated automatically**. Do not assume the old Chen launcher will
invoke this native implementation. Prepare an explicitly reviewed native
project configuration before using it. Checkpoint migration from that private
runtime or an old RC4 run is not certified by this addition; existing
checkpoint/provenance checks remain in force. Keep existing working projects
and their provenance intact during any separately planned migration.
