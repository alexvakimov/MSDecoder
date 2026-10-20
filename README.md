# MSDecoder: mass-spectrometry formula suggestions

MSDecoder is a Python command-line tool for ranking elemental-formula
candidates from a known precursor m/z and a positive-ion mass spectrum. It
combines precursor mass agreement, fragment subformula matches, common
neutral losses, a simple halogen isotope check, and a composition preference.
It produces candidate tables, reports, and spectrum plots. The current
implementation uses deterministic enumeration and heuristic scoring; it does
not train a machine-learning model or identify a unique molecular structure.

## Repository contents

```text
MSDecoder/
├── mass_spec_candidates.py       # Analysis functions and command-line interface
├── test_mass_spec_candidates.py # Automated regression and integration tests
├── requirements.txt             # Matplotlib for spectrum plots
├── requirements-structures.txt  # Base dependencies plus RDKit for structures
├── spectra/
│   └── slide02.txt … slide11.txt # Ten transcribed example peak lists
├── references/
│   └── slide_structures.csv     # Spectrum filename → known reference SMILES
├── results/
│   ├── summary.csv             # Saved batch ranking summary
│   └── slide02/ … slide11/      # Saved outputs for each example spectrum
├── README.md                   # Usage, scoring equations, and interpretation
└── LICENSE                     # GNU General Public License, version 3
```

| File or directory | Purpose |
| --- | --- |
| [mass_spec_candidates.py](mass_spec_candidates.py) | Reads spectra, enumerates bounded elemental formulae, ranks candidates, and writes outputs. Also supports a built-in demo and CSV or inline peaks. |
| [test_mass_spec_candidates.py](test_mass_spec_candidates.py) | Standard-library `unittest` suite covering formula enumeration, ranking, CLI settings, intensity weighting, and generated outputs. |
| [spectra/](spectra/) | Example input data. Each text file has metadata comments and `mz intensity` rows; metadata supplies the precursor, adduct, tolerance, and acquisition context. |
| [references/slide_structures.csv](references/slide_structures.csv) | Maps `spectrum_file` to `smiles` for optional drawings. These are supplied reference structures; they do not constrain or train the ranking algorithm. |
| [results/](results/) | Saved example analysis artifacts. Regenerate into a new output directory when comparing settings or code versions; saved results are not automatically refreshed by code changes. |
| [requirements.txt](requirements.txt) | Installs Matplotlib for PNG spectrum plots. Core formula ranking uses Python's standard library. |
| [requirements-structures.txt](requirements-structures.txt) | Includes the base requirements and adds RDKit for rendering supplied SMILES. |
| [LICENSE](LICENSE) | Contains the GNU GPL version 3 license text. |

The current repository does not include the original `Data.pptx` presentation.
The presentation is the recorded source of the transcriptions, not a runtime input: the script reads the supplied text files
and does not extract spectra or structures from PowerPoint.

### Example data and provenance

The files in `spectra/` are text peak lists transcribed from the spectra shown
on slides 2–11 of `Data.pptx`. Each file has comment-style acquisition metadata
followed by two columns, `mz intensity`. Peak m/z values are the labels visible
in the slide images; relative intensities are estimates from the plotted peak
heights. These are not the original instrument export files. Slides 5 and 7
are zoomed views of the spectra on slides 4 and 6, respectively.

### Analysis workflow and output files

1. Select a peak list and establish its precursor m/z and ion/adduct type.
2. Set mass tolerance and any independently justified element restrictions.
3. Run the script on one spectrum or a directory of text spectra.
4. Inspect the ranked formulae, mass errors, signal coverage, and spectrum plot.
5. Optionally supply reference SMILES to render known structures beside the results.

The text-file and batch modes write the following artifacts. Demo and
CSV/inline modes print rankings to the terminal instead.

| Output | Contents |
| --- | --- |
| `spectrum.png` | Stick plot normalized to the largest input intensity, with the selected precursor marked. |
| `candidates.csv` | Rank, formula, neutral exact mass, ppm error, composite score, fragment error, matched intensity fraction, and assigned fragment explanations. |
| `report.txt` | Input identity, analysis settings, candidate rankings, and reference SMILES when provided. |
| `structure.png` | Optional drawing of the supplied reference SMILES; it is not a structure predicted from the spectrum. |
| `summary.csv` | Batch-only overview with each spectrum's top candidate, diagnostics, alternative formulae, element settings, and optional reference SMILES. |

The saved example outputs are organized by input filename stem, such as
`spectra/slide02.txt` → `results/slide02/`. An existing output directory can
have files overwritten on a new run. Use a different `--output-dir` to retain
an earlier run. A high rank indicates agreement with the implemented scoring
rules, not confirmation of identity; the theory and limitations below explain
how to interpret each score.

## Run

Run the commands below from the repository root (the directory containing
`mass_spec_candidates.py`). Paths are relative, so the same commands work
after moving or cloning the repository. Use Python 3.10 or newer.

Install plotting support:

```sh
python -m pip install -r requirements.txt
```

Analyze one text spectrum:

```sh
python mass_spec_candidates.py --input spectra/slide02.txt --output-dir results/slide02
```

## Choosing elements

Formula enumeration has a configurable allowed-element whitelist. The default
is `C,H,N,O,S,P,Cl,Br`; the search only considers elements in that list. For
example, restrict a spectrum to carbon, hydrogen, nitrogen, oxygen, and
bromine, and require every candidate to contain bromine:

```sh
python mass_spec_candidates.py --input spectra/slide03.txt \
  --elements C,H,N,O,Br --require-elements Br \
  --output-dir results/slide03-brominated
```

`--elements` specifies what may occur, not what must occur. Use
`--require-elements` to require one or more selected elements in every result.
Required elements must also be in the allowed list. In text input files, set
the same options with metadata, for example:

```text
# elements: C,H,N,O,Br
# required_elements: Br
```

CLI options override the corresponding text-file metadata. At present the
supported element symbols are `C,H,N,O,S,P,Cl,Br`. Formula atom counts also
have built-in upper bounds; change `MAX_COUNTS` and `MAX_HYDROGENS` in the
script if a compound is outside those ranges.
Element restrictions also work with `--demo`, `--peaks`, and inline `--peak`
inputs.

Analyze all slide spectra and also render the separately recorded reference
structures from the presentation:

```sh
python -m pip install -r requirements-structures.txt
python mass_spec_candidates.py --batch spectra --structures references/slide_structures.csv --output-dir results
```

Each output directory contains `spectrum.png`, `candidates.csv`, and
`report.txt`. A batch run also writes `results/summary.csv`. When a reference
SMILES is provided, the spectrum's output directory also contains
`structure.png`. Reference structures are explicitly distinguished from
structures inferred by the spectrum algorithm.

For a new instrument export, convert or copy its peak list to a text file with
a metadata header and `mz intensity` columns. Provide `precursor_mz` and
`adduct` in the header; supported adduct settings are `proton`, `sodium`,
`ammonium`, and `none`. The precursor mass tolerance defaults to 25 ppm and can
be set with `tolerance_ppm` in the text header or `--tolerance-ppm` in any
input mode. An explicit CLI tolerance overrides the text header; otherwise
the header takes precedence over the 25 ppm default. The built-in `--demo` uses 50 ppm to match the rounded
two-decimal peak labels; the digitized text files likewise set their own wider
tolerances to account for rounding.

The precursor m/z is the selected molecular-ion peak measured by the mass
spectrometer, not a value calculated from the fragment peaks by this script.
For protonated chemical-ionization spectra, this is often the `[M+H]+` peak;
the input file's `precursor_mz` should be set to that observed peak and
`adduct: proton` tells the script to subtract the proton mass when estimating
the neutral molecule's mass. For an `[M]+` molecular ion, use `adduct: none`.
Choosing the correct precursor peak and ion/adduct type is an experimental
interpretation step. If the molecular ion is weak or absent, do not substitute
the base peak automatically: use another known molecular-ion/adduct peak, or
treat precursor-mass-based formula ranking as unavailable until a precursor
can be established.

## Tests

The [test_mass_spec_candidates.py](test_mass_spec_candidates.py) suite uses
Python's built-in `unittest` runner. Run it from the repository root after
installing the structure dependencies: the full suite generates both spectrum
plots and reference structure images, so it requires Matplotlib and RDKit.

```sh
python -m pip install -r requirements-structures.txt
python -m unittest -v test_mass_spec_candidates
```

The 14 tests cover:

- Formula enumeration across a wide hydrogen-count tolerance window and
  rejection of invalid masses or atom-count bounds.
- Allowed and required elements, including empty iterable overrides.
- CLI tolerance precedence, demo element restrictions, and conflicting inputs.
- The demo's top-ranked formula and the presence of each supplied spectrum's
  reference formula among its top ten candidates.
- The effect of peak intensity on fragment error and matched signal coverage.
- Single-file plots, reports, and CSV diagnostics, plus batch summary generation
  and reference structure rendering.

Output-producing tests use temporary directories, leaving the saved `results/`
artifacts untouched. The reference-formula tests check regression behavior on
the example data; they do not require every reference formula to rank first
or establish that the heuristic scores identify structures reliably.

## Theory: formula search and ranking

The implementation searches elemental formulae, not molecular graphs or
SMILES. The equations below describe the code as implemented. Higher scores
are better; scores are not probabilities. All ions are treated as singly
charged. Let $x_0$ be the observed precursor m/z, $a$ the selected adduct mass
offset, and $t$ the tolerance in ppm. Masses are expressed in Da.

### Neutral mass and formula enumeration

$$
m_0 = x_0-a, \qquad
m(\mathbf n)=\sum_{e} n_e\mu_e, \qquad
\delta_0=m_0\frac{t}{10^6}.
$$

Here $n_e$ is an element's atom count and $\mu_e$ its monoisotopic atomic
mass. A candidate must satisfy

$$
|m(\mathbf n)-m_0|\leq\delta_0, \qquad
\epsilon_{\mathrm{ppm}}=10^6\frac{m(\mathbf n)-m_0}{m_0}.
$$

The ppm denominator is the **neutral target mass**, not observed precursor
m/z. The code enumerates bounded heavy-element counts and all integer
hydrogen counts compatible with the remaining mass window. The allowed and
required element settings further restrict this search. Parent formulae must
contain carbon and pass the nonnegative integer double-bond-equivalent check:

$$
\mathrm{DBE}=\frac{2n_{\mathrm C}+2+n_{\mathrm N}+n_{\mathrm P}
-n_{\mathrm H}-n_{\mathrm{Cl}}-n_{\mathrm{Br}}}{2}
\in\{0,1,2,\ldots\}.
$$

One ring or double bond contributes one DBE; a triple bond contributes two.
This screening rule assumes trivalent N and P and ordinary closed-shell
neutral chemistry. It does not establish a valid structure. Fragment
subformulae are not subject to this rule, allowing radical compositions.
The `none` adduct uses zero offset; electron mass is neglected, which is an
approximation for a molecular radical cation at high mass accuracy.

### Composite score and mass-error terms

$$
S=-E_{\mathrm{prec}}-E_{\mathrm{frag}}+0.10L+0.20I+0.20P,
\qquad
E_{\mathrm{prec}}=\left(\frac{\epsilon_{\mathrm{ppm}}}{t}\right)^2.
$$

Let $\mathcal F$ contain the positive-intensity peaks remaining after the
precursor isotope mask. For peak $i$, write observed m/z as $x_i$ and signal
intensity as $J_i$. The mask excludes peaks near offsets 0, 1.003355, and
1.99705 Da from the precursor, using a window of
$\max(x_0t/10^6,0.01)$ Da. It is a positional heuristic, not full deisotoping:
fragment isotopes and higher precursor isotopes can remain.

$$
w_i=\frac{J_i}{\sum_{j\in\mathcal F}J_j},\qquad
\delta_i=x_i\frac{t}{10^6},\qquad
\widehat{x}_i=m(\mathbf f_i)+a.
$$

The fragment composition $\mathbf f_i$ must fit componentwise within the
parent formula. Using the precursor's offset $a$ for every fragment assumes
that the ion/adduct type is retained. Adduct exchange and multiple charges
are not modeled.

$$
q_i=
\begin{cases}
\min\!\left(1,\left(\dfrac{\widehat{x}_i-x_i}{\delta_i}\right)^2\right),
&\text{if assigned a subformula within tolerance},\\
1,&\text{if unmatched},
\end{cases}
\qquad
E_{\mathrm{frag}}=\sum_{i\in\mathcal F}w_iq_i.
$$

Each peak and each unique subformula can be assigned at most once. The parent
formula itself is excluded as a fragment assignment. The greedy matcher
selects available assignments in descending order of improvement
$w_i(1-q_i)$; it does not guarantee a globally optimal assignment.

The reported signal coverage is

$$
C_{\mathrm{matched}}=\sum_{i\in\mathcal A}w_i,
$$

where $\mathcal A$ is the assigned subset. Both coverage and fragment error
lie between zero and one. If no eligible peaks exist, the code sets
$E_{\mathrm{frag}}=1$ and $C_{\mathrm{matched}}=0$ by convention; that default
is not evidence of a failed fragmentation mechanism. Peaks above the
precursor remain in $\mathcal F$ unless masked and generally receive unmatched
penalties. In these full-scan inputs, such peaks may belong to other species.
Thus coverage describes the supplied eligible signal, not necessarily the
fraction of signal originating from the selected precursor.

### Neutral-loss support, isotope fit, and composition preference

For each eligible peak, $b_i=1$ if at least one listed neutral loss fits the
parent's atom counts and the observed precursor-to-peak difference; otherwise
$b_i=0$. Each peak contributes only once:

$$
L=\sum_{i\in\mathcal F}w_i b_i,\qquad
\left|(x_0-x_i)-m(\boldsymbol\ell)\right|
\leq\max\!\left(0.02,\,m(\boldsymbol\ell)\frac{t}{10^6}\right).
$$

The supported loss compositions are $\mathrm{CH_4O}$, $\mathrm{H_2O}$,
$\mathrm{CO}$, $\mathrm{CO_2}$, $\mathrm{NH_3}$, and $\mathrm{HCN}$.
The implementation uses their tabulated, rounded exact masses. This check
assumes identical charge/adduct offsets cancel and does not require the loss
to agree with the selected fragment assignment. It is weak relational
evidence, not a fragmentation tree. With no eligible peaks, $L=0$.

For a halogen-containing candidate with usable precursor and M+2 signals:

$$
r_{\mathrm{obs}}=\frac{J_{M+2}}{J_M},\qquad
I=\max\!\left(0,1-\frac{|r_{\mathrm{obs}}-r_{\mathrm{expected}}|}{0.5}\right),
\qquad
r_{\mathrm{expected}}=
\begin{cases}
0.97,&n_{\mathrm{Br}}>0,\\
0.32,&n_{\mathrm{Br}}=0\text{ and }n_{\mathrm{Cl}}>0.
\end{cases}
$$

The strongest peaks within $\max(0.02,x_0t/10^6)$ Da of $x_0$ and
$x_0+1.99705$ are used. Otherwise $I=0$, including candidates without halogens.
This is a **single-halogen approximation**: it ignores halogen multiplicity,
mixed Cl/Br compositions, and other M+2 contributions. It does not compare a
complete theoretical isotope envelope.

Let $h$ be the total N, O, S, P, Cl, and Br count. The implemented composition
preference is

$$
P=\frac{0.5\min(n_{\mathrm C},20)+0.75\min(\mathrm{DBE},12)
-2\max(0,h-n_{\mathrm C}-1)}{20}.
$$

This is an uncalibrated preference for carbon and moderate unsaturation; it
can be negative and is not a probability. Score coefficients are heuristic.
Fragment matching and neutral-loss support can reward the same peak, so these
terms should not be interpreted as independent evidence or log likelihoods.

## What does intensity weighting mean?

**Intensity weighting is meaningful as an emphasis on explaining abundant
observed ion signal. It is not, by itself, a statistical model of mass error
or fragmentation yield.** Intensity reflects detected ions, affected by
formation pathways, subsequent fragmentation, ion survival, transmission,
and detector response. Neutral products are not observed as ordinary ion
peaks. Relative peak heights therefore need not equal the fractions of
precursor molecules choosing particular primary cleavage pathways.

A schematic response model, **not implemented here**, is

$$
J_i\approx R_i N_i^{\mathrm{surviving}},
$$

where $R_i$ represents instrumental response and
$N_i^{\mathrm{surviving}}$ is the number of ions reaching detection through all
relevant pathways. The script has neither a prediction for those populations
nor calibrated response factors. In this repository, intensities are also
visual estimates from plots, not raw ion counts.

For example, an exact match with intensity 90 and an unmatched peak with
intensity 10 gives $E_{\mathrm{frag}}=0.10$ and coverage 90%. Reversing the
intensities gives error 0.90 and coverage 10%. The program is choosing to
prioritize abundant signal; it has not demonstrated that the first peak's
mass is nine times more precise or that its assignment is chemically correct.
Multiplying all intensities by a common factor leaves these results unchanged.

Statistical weighted least squares would instead need a justified mass-error
variance model, with weights proportional to $1/\sigma_i^2$. A ppm tolerance
is an acceptance window, not automatically a standard deviation. Greater
signal can improve centroid precision, but the relationship requires
instrument-specific calibration; systematic calibration error and rounded
labels remain. Linear intensity weighting is especially uncertain for the
manually transcribed data here.

Possible future comparisons include uniform weights or square-root intensity
weights:

$$
w_i^{(\alpha)}=\frac{J_i^\alpha}{\sum_{j\in\mathcal F}J_j^\alpha},
\qquad \alpha\in\{0,\tfrac12,1\}.
$$

Only $\alpha=1$ is implemented. Choosing an alternative requires validation
on independent spectra rather than tuning to make the supplied reference
formulae rank first. Predicting fragment intensities would be a separate
model, usually requiring structure, ionization/activation conditions, and
training data or a mechanistic calculation.

## Relationships between peaks and mass conservation

### Why a complementary mass need not appear as a peak

A common singly charged dissociation is

$$
\mathrm{P}^{+}\longrightarrow\mathrm{F}^{+}+\mathrm{L}^{0}.
$$

Atoms and charge are conserved, but the complementary product is neutral.
For this channel, with both ions singly charged, their m/z difference gives
the neutral-loss mass to the usual chemical mass accuracy:

$$
m(\mathrm{L}^{0})\approx m(\mathrm{P}^{+})-m(\mathrm{F}^{+})=x_0-x_F.
$$

Thus a precursor at m/z 200 yielding a product at m/z 182 can indicate an
approximately 18 Da neutral loss; it does **not** require a peak at m/z 18.
The neutral molecular mass $M$ also differs from the precursor ion mass when
an adduct is present. Subtracting a raw fragment m/z from $M$ without accounting
for charge and adduct is generally incorrect.

For a channel yielding two charged products, charge as well as composition
must balance. For example, a doubly charged precursor can produce two singly
charged products. With consistent ion masses and no additional products:

$$
z_P=z_F+z_G,\qquad z_Px_P\approx z_Fx_F+z_Gx_G.
$$

This is not the singly charged model used by the script. Peaks representing
alternative charge-retention channels can also arise from different precursor
ions in the ensemble; their coexistence does not imply two positive fragments
came from one singly charged dissociation or that their intensities are equal.

### Sequential losses and fragmentation trees

Peaks can be chemically related. For example,
$\mathrm{P}^{+}\to\mathrm{F}^{+}+\mathrm{L}_1^0$ followed by
$\mathrm{F}^{+}\to\mathrm{G}^{+}+\mathrm{L}_2^0$ gives a sequence of losses.
A possible extension would assign formulae to nodes and neutral losses to
edges of a directed fragmentation tree. For an edge $u\to v$, under a shared
singly charged ion/adduct convention, require

$$
\boldsymbol\ell_{uv}=\mathbf f_u-\mathbf f_v\geq\mathbf 0,
\qquad \boldsymbol\ell_{uv}\ne\mathbf 0,\qquad
\left|(x_u-x_v)-m(\boldsymbol\ell_{uv})\right|\leq\delta_{uv}.
$$

Formula subtraction enforces atom balance; charge balance and plausible loss
chemistry are additional requirements. The edge tolerance must account for
both measurements. For calibrated independent uncertainties,
$\sigma_{uv}^2=\sigma_u^2+\sigma_v^2$; more generally,

$$
\sigma_{uv}^2=\sigma_u^2+\sigma_v^2-2\operatorname{Cov}(x_u,x_v).
$$

The code's current neutral-loss window is a heuristic, not this uncertainty
model. Rearrangements can transfer hydrogen, and changes in adduct or charge
require explicit ion-composition bookkeeping before applying these equations.
Even balanced formula differences do not prove a bond-breaking route.

A future joint score could combine mass fit with a validated tree score,
$S_{\mathrm{joint}}=S+\lambda T$, but $T$ and $\lambda$ are **not implemented**.
It should select a consistent subset of edges rather than reward every
possible pair, avoid double-counting existing loss evidence, and allow missing
intermediates. For each individual cleavage atom counts balance, but one must
not sum all observed fragment formulae and require that sum to fit one parent:
the spectrum aggregates many precursor ions and alternative pathways.

Fragmentation trees are an established way to encode these relationships;
see the primary study [Fragmentation trees reloaded](https://pmc.ncbi.nlm.nih.gov/articles/PMC4736045/)
and the [SIRIUS methods documentation](https://v6.docs.sirius-ms.io/methods-background/).
Those methods use more developed models than this script's independent
subformula matching and precursor-only common-loss bonus.

### Can fragments combine, and can correlations be imposed?

Successive **losses** can add: the total precursor-to-final-product mass
loss is the sum along a pathway. This does not mean two observed fragment
ions recombine. Ion–molecule association, clustering, or adduct formation can
produce heavier ions under suitable conditions, but a numerical sum of two
peak masses is insufficient evidence for such a reaction. Any association
hypothesis needs neutral/ionic reactants, charge and atom balance, and
experimental support; it is not a general fragmentation constraint.

Isotope envelopes, related neutral-loss series, and shared precursor origin
can provide useful joint constraints. Statistical correlation is different
from a chemically plausible mass relation: a single peak list cannot estimate
intensity covariance across observations. Replicate spectra, chromatographic
profiles, collision-energy series, or MS/MS and MS3 experiments can provide
additional evidence, though co-elution or co-isolation can still confound
shared-origin inference. Equal intensities or fixed complementary-peak ratios
should not be imposed without a specific validated mechanism and response
model.

For these files, limited scan ranges, rounded labels, estimated intensities,
and the absence of isolated-precursor MS/MS evidence make strong pathway
constraints premature. Slides 5 and 7 are zoom views, not independent
replicates of slides 4 and 6. The current implementation imposes parent
subformula bounds, unique assignments, a precursor-loss bonus, and a rough
halogen isotope check; it does **not** fit a fragmentation tree, fragment
recombination, intensity covariance, or predicted fragment yields.

For ion, neutral-loss, and intensity terminology, see the
[IUPAC recommendations for mass spectrometry](https://publications.iupac.org/pac/85/7/1515/index.html).

## Interpretation limits

Rounded two-decimal labels make formula matches ambiguous. Narrowing the
allowed elements using independent sample knowledge can help, but an incorrect
whitelist excludes the true formula. Reference SMILES and drawings are
supplied annotations, not inferred structures or inputs to candidate scoring.
The reference mapping in `references/slide_structures.csv` records no Br for
slide 2 and a brominated analogue for slide 3; this documentation does not
require reading `Data.pptx`.

The tool does not search a structure database or identify an isomer. Better
mass calibration, raw spectra, isolated-precursor MS/MS, and a validated
fragmentation model would be needed for stronger structural conclusions.
