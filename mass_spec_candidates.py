#!/usr/bin/env python3
"""Rank elemental-formula candidates from a precursor and positive-ion MS peaks.

This is a formula-ranking aid, not an isomer-identification engine. Without a
structure database or MS/MS interpretation, spectra generally cannot determine
a unique molecular structure.

Demo on the first spectrum in Data.pptx:
    python mass_spec_candidates.py --demo

Analyze a spectrum from CSV (header columns: mz,intensity):
    python mass_spec_candidates.py --precursor 257.16 --peaks peaks.csv

Restrict formula candidates to carbon, hydrogen, nitrogen, oxygen, and bromine,
and require bromine in every candidate:
    python mass_spec_candidates.py --input spectra/slide03.txt \\
        --elements C,H,N,O,Br --require-elements Br

Or pass peaks inline:
    python mass_spec_candidates.py --precursor 257.16 \\
        --peak 225.14:85 --peak 257.16:100
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping


PROTON_MASS = 1.007276466621
ADDUCT_MASSES = {
    "proton": PROTON_MASS,
    "sodium": 22.989218,
    "ammonium": 18.033823,
    "none": 0.0,
}
ATOMIC_MASSES = {
    "C": 12.0,
    "H": 1.00782503223,
    "N": 14.00307400443,
    "O": 15.99491461957,
    "S": 31.9720711744,
    "P": 30.97376199842,
    "Cl": 34.968852682,
    "Br": 78.9183376,
}
ELEMENT_ORDER = ("C", "H", "N", "O", "S", "P", "Cl", "Br")
HEAVY_ELEMENTS = ("C", "N", "O", "S", "P", "Cl", "Br")
DEFAULT_ELEMENTS = ELEMENT_ORDER
MAX_COUNTS = {"C": 30, "N": 10, "O": 12, "S": 4, "P": 4, "Cl": 4, "Br": 2}
MAX_HYDROGENS = 60

# Values are neutral-loss formulae and exact masses. Their presence can support
# a candidate composition, but is not proof of a particular bond or structure.
COMMON_NEUTRAL_LOSSES = (
    ({"C": 1, "H": 4, "O": 1}, 32.026215),
    ({"H": 2, "O": 1}, 18.010565),
    ({"C": 1, "O": 1}, 27.994915),
    ({"C": 1, "O": 2}, 43.989830),
    ({"N": 1, "H": 3}, 17.026549),
    ({"C": 1, "H": 1, "N": 1}, 27.010899),
)


@dataclass(frozen=True)
class Peak:
    """Observed m/z and nonnegative detector signal in arbitrary relative units.

    Intensity weights observed signal coverage; it is not a calibrated
    fragmentation branching probability or a mass-error precision estimate.
    """

    mz: float
    intensity: float


@dataclass(frozen=True)
class Candidate:
    """Formula with a heuristic score (larger is better) and fit diagnostics.

    ``exact_mass`` is neutral mass in Da; ``ppm_error`` uses the neutral target
    as denominator. ``fragment_mse`` is tolerance-normalized weighted error,
    not an error in Da squared. ``matched_intensity_fraction`` measures signal
    coverage, not the fraction of precursor molecules following a pathway.
    """

    formula: str
    exact_mass: float
    ppm_error: float
    score: float
    explained_peaks: tuple[str, ...]
    fragment_mse: float = 1.0
    matched_intensity_fraction: float = 0.0


def formula_mass(counts: Mapping[str, int]) -> float:
    """Calculate a neutral monoisotopic mass from elemental atom counts."""
    return sum(ATOMIC_MASSES[element] * count for element, count in counts.items())


def format_formula(counts: Mapping[str, int]) -> str:
    """Format atom counts in the script's conventional Hill-like element order."""
    return "".join(
        element + (str(counts[element]) if counts[element] != 1 else "")
        for element in ELEMENT_ORDER
        if counts.get(element, 0) > 0
    )


def _is_plausible_formula(counts: Mapping[str, int]) -> bool:
    """Require carbon and nonnegative integer (2C + 2 + N + P - H - Cl - Br)/2.

    This neutral closed-shell screening rule assumes trivalent N and P; it
    neither establishes a molecular graph nor applies to every ion chemistry.
    """
    carbon = counts.get("C", 0)
    hydrogen = counts.get("H", 0)
    nitrogen = counts.get("N", 0)
    phosphorus = counts.get("P", 0)
    halogens = counts.get("Cl", 0) + counts.get("Br", 0)
    dbe_twice = 2 * carbon + 2 + nitrogen + phosphorus - hydrogen - halogens
    return carbon > 0 and dbe_twice >= 0 and dbe_twice % 2 == 0


def _composition_prior(counts: Mapping[str, int]) -> float:
    """Return an uncalibrated preference for carbon and moderate unsaturation.

    The value is 0.5*min(C,20) + 0.75*min(DBE,12) minus twice the excess
    heteroatom count above C+1. It is not a learned log prior or probability.
    """
    carbon = counts.get("C", 0)
    heteroatoms = sum(counts.get(element, 0) for element in ("N", "O", "S", "P", "Cl", "Br"))
    dbe_twice = (
        2 * carbon
        + 2
        + counts.get("N", 0)
        + counts.get("P", 0)
        - counts.get("H", 0)
        - counts.get("Cl", 0)
        - counts.get("Br", 0)
    )
    dbe = max(0, dbe_twice / 2)
    hetero_excess = max(0, heteroatoms - carbon - 1)
    return 0.5 * min(carbon, 20) + 0.75 * min(dbe, 12) - 2.0 * hetero_excess


def enumerate_formulas(
    target_mass: float,
    tolerance_da: float,
    max_counts: Mapping[str, int] = MAX_COUNTS,
    hydrogen_limit: int = MAX_HYDROGENS,
    allowed_elements: Iterable[str] = DEFAULT_ELEMENTS,
    check_plausibility: bool = True,
) -> list[tuple[dict[str, int], float]]:
    """Enumerate neutral formulae within a mass tolerance and element whitelist.

    All integer ``H`` counts within the residual mass window are tried; if H is absent, only
    hydrogen-free formulae are considered. Heavy-element maxima remain bounded
    by ``max_counts`` to keep the search finite. Set ``check_plausibility`` to
    false when enumerating fragment-ion compositions that may be radicals.
    """
    if (
        not math.isfinite(target_mass)
        or not math.isfinite(tolerance_da)
        or target_mass <= 0
        or tolerance_da < 0
    ):
        raise ValueError("target mass must be finite and positive and tolerance finite and nonnegative")
    if any(
        not isinstance(bound, int) or isinstance(bound, bool) or bound < 0
        for bound in (*max_counts.values(), hydrogen_limit)
    ):
        raise ValueError("atom-count bounds must be nonnegative integers")

    element_set = set(parse_elements(allowed_elements))
    heavy_elements = tuple(element for element in HEAVY_ELEMENTS if element in element_set)
    use_hydrogen = "H" in element_set
    matches: list[tuple[dict[str, int], float]] = []

    def walk(index: int, partial_mass: float, counts: dict[str, int]) -> None:
        """Recursively try bounded counts for each permitted heavy element."""
        if index == len(heavy_elements):
            residual = target_mass - partial_mass
            # A wide tolerance can admit more than the nearest hydrogen count.
            # Expand the integer interval slightly, then apply the exact mass
            # check, so floating-point rounding cannot drop a boundary count.
            if use_hydrogen:
                lower = max(0, math.floor((residual - tolerance_da) / ATOMIC_MASSES["H"]))
                upper = min(hydrogen_limit, math.ceil((residual + tolerance_da) / ATOMIC_MASSES["H"]))
                hydrogen_counts = range(lower, upper + 1)
            else:
                hydrogen_counts = (0,)
            for h_count in hydrogen_counts:
                formula_counts = {**counts, "H": h_count}
                exact_mass = formula_mass(formula_counts)
                if (
                    abs(exact_mass - target_mass) <= tolerance_da
                    and (not check_plausibility or _is_plausible_formula(formula_counts))
                ):
                    matches.append((formula_counts, exact_mass))
            return

        element = heavy_elements[index]
        max_count = max_counts.get(element, 0)
        atomic_mass = ATOMIC_MASSES[element]
        for count in range(max_count + 1):
            next_mass = partial_mass + count * atomic_mass
            if next_mass > target_mass + tolerance_da:
                break
            counts[element] = count
            walk(index + 1, next_mass, counts)
        counts.pop(element, None)

    walk(0, 0.0, {})
    return matches


def _within_counts(smaller: Mapping[str, int], larger: Mapping[str, int]) -> bool:
    """Return whether every atom in ``smaller`` is available in ``larger``."""
    return all(count <= larger.get(element, 0) for element, count in smaller.items())


def _fragment_formula_matches(
    parent: Mapping[str, int],
    fragment_mass: float,
    tolerance_da: float,
    ion_mass_offset: float,
    allowed_elements: Iterable[str] = DEFAULT_ELEMENTS,
) -> list[tuple[dict[str, int], float]]:
    """Return mass-feasible subformulae and neutral atomic-sum masses in Da.

    Subtract ``ion_mass_offset`` from the observed singly charged fragment
    m/z, using the same offset as the precursor. This assumes adduct retention;
    it does not model adduct exchange, bond connectivity, or fragment yields.
    Atom counts cannot exceed the parent, but neutral closed-shell DBE is not
    enforced. A match alone does not establish a fragmentation pathway.
    """
    target_neutral_mass = fragment_mass - ion_mass_offset
    if target_neutral_mass <= 0:
        return []
    allowed_set = set(allowed_elements)
    bounds = {
        element: parent.get(element, 0)
        for element in HEAVY_ELEMENTS
        if element in allowed_set
    }
    return [
        (counts, exact_mass)
        for counts, exact_mass in enumerate_formulas(
            target_neutral_mass,
            tolerance_da,
            bounds,
            parent.get("H", 0),
            allowed_set,
            check_plausibility=False,
        )
        if _within_counts(counts, parent)
    ]


def _loss_supported(parent: Mapping[str, int], loss: Mapping[str, int]) -> bool:
    """Check whether a proposed neutral loss fits within a parent composition."""
    return _within_counts(loss, parent)


def _is_precursor_isotope_peak(
    peak_mz: float, precursor_mz: float, tolerance_ppm: float
) -> bool:
    """Mask peaks near precursor offsets 0, 1.003355, and 1.99705 Da.

    The window is max(precursor_mz*tolerance_ppm/1e6, 0.01 Da).
    This positional heuristic does not establish isotope identity or remove
    fragment isotope envelopes and higher precursor isotopes.
    """
    tolerance = max(
        precursor_mz * tolerance_ppm / 1_000_000,
        0.01,
    )
    return any(
        abs(peak_mz - (precursor_mz + isotope_shift)) <= tolerance
        for isotope_shift in (0.0, 1.003355, 1.99705)
    )


def _score_fragment_spectrum(
    parent: Mapping[str, int],
    peaks: Iterable[Peak],
    precursor_mz: float,
    ion_mass_offset: float,
    tolerance_ppm: float,
    allowed_elements: Iterable[str],
) -> tuple[float, float, tuple[str, ...]]:
    """Return intensity-weighted normalized fragment MSE, coverage, and matches.

    Each experimental peak can be assigned to at most one unique elemental
    sub-formula. Peaks with no match inside the mass tolerance receive the
    maximum normalized squared error of 1.0. Precursor isotope peaks are
    excluded because they are scored separately.

    For eligible positive-intensity peaks, w_i = intensity_i/sum(intensity)
    and matched error is ((predicted_mz-observed_mz)/tolerance_Da)^2.
    Greedy assignment maximizes the next weighted improvement over error 1;
    it need not find the globally best assignment. Return sum(w_i*error_i),
    assigned intensity fraction, and textual matches. With no eligible peaks,
    return (1, 0, ()).

    Intensity emphasizes abundant detected signal, not inverse measurement
    variance or predicted fragmentation yield. Above-precursor peaks remain
    in the denominator and generally receive unmatched penalties. No
    fragment-to-fragment losses, complementary pairs, or covariance are fit.
    """
    fragment_peaks = [
        peak
        for peak in peaks
        if peak.intensity > 0
        and not _is_precursor_isotope_peak(peak.mz, precursor_mz, tolerance_ppm)
    ]
    total_intensity = sum(peak.intensity for peak in fragment_peaks)
    if total_intensity <= 0:
        return 1.0, 0.0, ()

    peak_weights = [peak.intensity / total_intensity for peak in fragment_peaks]
    errors = [1.0] * len(fragment_peaks)
    assignments: list[tuple[float, float, int, str, float]] = []
    for peak_index, (peak, weight) in enumerate(zip(fragment_peaks, peak_weights)):
        tolerance_da = peak.mz * tolerance_ppm / 1_000_000
        if tolerance_da <= 0:
            continue
        matches = _fragment_formula_matches(
            parent,
            peak.mz,
            tolerance_da,
            ion_mass_offset,
            allowed_elements,
        )
        for fragment_counts, fragment_mass in matches:
            if format_formula(fragment_counts) == format_formula(parent):
                continue
            observed_neutral_mass = peak.mz - ion_mass_offset
            normalized_error = (fragment_mass - observed_neutral_mass) / tolerance_da
            squared_error = normalized_error * normalized_error
            benefit = weight * (1.0 - squared_error)
            assignments.append(
                (
                    -benefit,
                    squared_error,
                    peak_index,
                    format_formula(fragment_counts),
                    normalized_error,
                )
            )

    assigned_peaks: set[int] = set()
    assigned_fragments: set[str] = set()
    explanations: list[str] = []
    for _, squared_error, peak_index, fragment_formula, normalized_error in sorted(assignments):
        if peak_index in assigned_peaks or fragment_formula in assigned_fragments:
            continue
        assigned_peaks.add(peak_index)
        assigned_fragments.add(fragment_formula)
        errors[peak_index] = min(1.0, squared_error)
        peak = fragment_peaks[peak_index]
        explanations.append(
            f"{peak.mz:.2f}:{fragment_formula}"
            f"({normalized_error:+.2f} tolerance units)"
        )

    fragment_mse = sum(
        weight * error for weight, error in zip(peak_weights, errors)
    )
    matched_intensity_fraction = sum(
        peak_weights[index] for index in assigned_peaks
    )
    return fragment_mse, matched_intensity_fraction, tuple(explanations)


def _neutral_loss_support(
    parent: Mapping[str, int],
    peaks: Iterable[Peak],
    precursor_mz: float,
    tolerance_ppm: float,
) -> float:
    """Return eligible signal fraction supporting a listed precursor loss.

    Require loss atoms to fit the parent and precursor_mz-peak.mz to agree
    within max(0.02 Da, loss_mass*tolerance_ppm/1e6). Count each peak once.
    Equal charge/adduct offsets are assumed to cancel. This is a separate
    heuristic from fragment assignment, not a joint pathway or yield model;
    it may reward the same evidence already used in fragment scoring.
    """
    fragment_peaks = [
        peak
        for peak in peaks
        if peak.intensity > 0
        and not _is_precursor_isotope_peak(peak.mz, precursor_mz, tolerance_ppm)
    ]
    total_intensity = sum(peak.intensity for peak in fragment_peaks)
    if total_intensity <= 0:
        return 0.0

    supported_intensity = 0.0
    for peak in fragment_peaks:
        loss_mass = precursor_mz - peak.mz
        for loss_formula, exact_loss_mass in COMMON_NEUTRAL_LOSSES:
            loss_tolerance = max(0.02, exact_loss_mass * tolerance_ppm / 1_000_000)
            if (
                _loss_supported(parent, loss_formula)
                and abs(loss_mass - exact_loss_mass) <= loss_tolerance
            ):
                supported_intensity += peak.intensity
                break
    return supported_intensity / total_intensity


def _halogen_isotope_fit(
    counts: Mapping[str, int],
    peaks: Iterable[Peak],
    precursor_mz: float,
    tolerance_ppm: float,
) -> float:
    """Return max(0, 1-abs(observed_ratio-expected_ratio)/0.5).

    Use the strongest peaks near precursor and precursor+1.99705 Da, with
    max(0.02 Da, precursor_mz*tolerance_ppm/1e6) windows. Expected ratio is
    0.97 if any Br occurs, otherwise 0.32 for Cl. Return zero without halogens
    or usable peaks. This single-halogen approximation ignores atom
    multiplicity, mixed halogens, and other isotope contributions.
    """
    if not (counts.get("Cl", 0) or counts.get("Br", 0)):
        return 0.0

    isotope_tolerance = max(0.02, precursor_mz * tolerance_ppm / 1_000_000)
    peak_list = tuple(peaks)
    precursor_peak = max(
        (
            peak
            for peak in peak_list
            if abs(peak.mz - precursor_mz) <= isotope_tolerance
        ),
        key=lambda peak: peak.intensity,
        default=None,
    )
    isotope_peak = max(
        (
            peak
            for peak in peak_list
            if abs(peak.mz - precursor_mz - 1.99705) <= isotope_tolerance
        ),
        key=lambda peak: peak.intensity,
        default=None,
    )
    if precursor_peak is None or isotope_peak is None or precursor_peak.intensity <= 0:
        return 0.0

    observed_ratio = isotope_peak.intensity / precursor_peak.intensity
    expected_ratio = 0.97 if counts.get("Br", 0) else 0.32
    return max(0.0, 1.0 - abs(observed_ratio - expected_ratio) / 0.5)


def parse_elements(value: str | Iterable[str]) -> tuple[str, ...]:
    """Validate and normalize element symbols from a comma-separated list."""
    raw_elements = value.split(",") if isinstance(value, str) else value
    elements = tuple(dict.fromkeys(element.strip() for element in raw_elements if element.strip()))
    unknown = sorted(set(elements) - set(ATOMIC_MASSES))
    if unknown:
        raise ValueError(
            f"unsupported element(s): {', '.join(unknown)}; choose from "
            f"{', '.join(ELEMENT_ORDER)}"
        )
    if not elements:
        raise ValueError("at least one allowed element must be specified")
    return tuple(element for element in ELEMENT_ORDER if element in elements)


def rank_candidates(
    precursor_mz: float,
    peaks: Iterable[Peak],
    adduct: str = "proton",
    tolerance_ppm: float = 25.0,
    limit: int = 10,
    allowed_elements: Iterable[str] = DEFAULT_ELEMENTS,
    required_elements: Iterable[str] = (),
) -> list[Candidate]:
    """Rank precursor formulae by mass fit, fragment evidence, and composition.

    ``allowed_elements`` constrains the formula search; it does not assert that
    every listed element must occur. Hydrogen is included by default and can be
    omitted to restrict candidates to hydrogen-free compositions.
    ``required_elements`` can further require at least one atom of selected
    elements in each parent formula.

    Score = -precursor_mse - fragment_mse + 0.10*loss_support
    + 0.20*isotope_fit + 0.20*composition_prior/20. Precursor MSE is squared
    ppm error divided by squared ppm tolerance, relative to neutral mass.
    Scores are heuristic, not probabilities or independent log likelihoods.
    All ions are treated as singly charged with the same adduct offset;
    ``none`` uses zero offset and neglects electron mass. No structure,
    predicted intensity, fragmentation tree, or peak covariance is inferred.
    """
    elements = parse_elements(allowed_elements)
    required_values = (
        required_elements if isinstance(required_elements, str) else tuple(required_elements)
    )
    required = parse_elements(required_values) if required_values else ()
    if not set(required) <= set(elements):
        raise ValueError("required elements must also be included in allowed elements")
    if adduct not in ADDUCT_MASSES:
        raise ValueError(f"unsupported adduct {adduct!r}; choose from {', '.join(ADDUCT_MASSES)}")
    if (
        not math.isfinite(precursor_mz)
        or not math.isfinite(tolerance_ppm)
        or precursor_mz <= ADDUCT_MASSES[adduct]
        or tolerance_ppm <= 0
    ):
        raise ValueError("precursor m/z must exceed the adduct mass and tolerance must be positive")
    if limit < 1:
        raise ValueError("limit must be at least 1")

    peak_list = tuple(peaks)
    if any(
        not math.isfinite(peak.mz)
        or not math.isfinite(peak.intensity)
        or peak.mz <= 0
        or peak.intensity < 0
        for peak in peak_list
    ):
        raise ValueError("peak m/z must be finite and positive; intensity must be finite and nonnegative")
    neutral_mass = precursor_mz - ADDUCT_MASSES[adduct]
    tolerance_da = neutral_mass * tolerance_ppm / 1_000_000
    formula_matches = [
        (counts, exact_mass)
        for counts, exact_mass in enumerate_formulas(
            neutral_mass, tolerance_da, allowed_elements=elements
        )
        if all(counts.get(element, 0) > 0 for element in required)
    ]
    if not formula_matches:
        return []

    ranked: list[Candidate] = []
    for counts, exact_mass in formula_matches:
        ppm_error = (exact_mass - neutral_mass) / neutral_mass * 1_000_000
        precursor_mse = (ppm_error / tolerance_ppm) ** 2
        fragment_mse, matched_intensity, explained = _score_fragment_spectrum(
            counts,
            peak_list,
            precursor_mz,
            ADDUCT_MASSES[adduct],
            tolerance_ppm,
            elements,
        )
        loss_support = _neutral_loss_support(
            counts, peak_list, precursor_mz, tolerance_ppm
        )
        isotope_fit = _halogen_isotope_fit(
            counts, peak_list, precursor_mz, tolerance_ppm
        )
        composition_prior = _composition_prior(counts) / 20.0
        score = (
            -precursor_mse
            - fragment_mse
            + 0.10 * loss_support
            + 0.20 * isotope_fit
            + 0.20 * composition_prior
        )

        ranked.append(
            Candidate(
                format_formula(counts),
                exact_mass,
                ppm_error,
                score,
                explained,
                fragment_mse,
                matched_intensity,
            )
        )

    ranked.sort(key=lambda candidate: (-candidate.score, abs(candidate.ppm_error), candidate.formula))
    return ranked[:limit]


def _parse_peak(value: str) -> Peak:
    """Parse one command-line peak supplied as ``m/z:intensity``."""
    try:
        mz_text, intensity_text = value.split(":", maxsplit=1)
        peak = Peak(float(mz_text), float(intensity_text))
    except (ValueError, TypeError) as error:
        raise argparse.ArgumentTypeError("peak must be specified as m/z:intensity") from error
    if not math.isfinite(peak.mz) or not math.isfinite(peak.intensity) or peak.mz <= 0 or peak.intensity < 0:
        raise argparse.ArgumentTypeError("peak m/z must be positive and intensity nonnegative")
    return peak


def _read_csv(path: Path) -> list[Peak]:
    """Read a CSV peak list containing ``mz`` and ``intensity`` columns."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not {"mz", "intensity"} <= set(reader.fieldnames):
            raise ValueError("peak CSV must have 'mz' and 'intensity' columns")
        peaks = []
        for line_number, row in enumerate(reader, start=2):
            try:
                peaks.append(Peak(float(row["mz"]), float(row["intensity"])))
            except (TypeError, ValueError) as error:
                raise ValueError(f"invalid peak values on CSV line {line_number}") from error
    return peaks


def read_text_spectrum(path: Path) -> tuple[dict[str, str], list[Peak]]:
    """Read comment metadata and whitespace-, tab-, or comma-separated peaks."""
    metadata: dict[str, str] = {}
    peaks: list[Peak] = []
    in_data = False
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                entry = line[1:].strip()
                if ":" in entry:
                    key, value = entry.split(":", maxsplit=1)
                    metadata[key.strip().lower()] = value.strip()
                continue
            columns = line.replace(",", " ").split()
            if not in_data and [column.lower() for column in columns[:2]] == ["mz", "intensity"]:
                in_data = True
                continue
            in_data = True
            if len(columns) < 2:
                raise ValueError(f"{path}:{line_number}: expected m/z and intensity columns")
            try:
                peak = Peak(float(columns[0]), float(columns[1]))
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid numeric peak") from error
            if not math.isfinite(peak.mz) or not math.isfinite(peak.intensity) or peak.mz <= 0 or peak.intensity < 0:
                raise ValueError(f"{path}:{line_number}: m/z must be positive and intensity nonnegative")
            peaks.append(peak)
    if not peaks:
        raise ValueError(f"{path}: no peak rows found")
    return metadata, peaks


def _draw_spectrum(path: Path, title: str, precursor_mz: float, peaks: Iterable[Peak]) -> None:
    """Save a relative-intensity stick plot and precursor marker as a PNG."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError("plotting requires matplotlib; install requirements.txt") from error

    peak_list = tuple(peaks)
    maximum = max(peak.intensity for peak in peak_list)
    figure, axis = plt.subplots(figsize=(11, 5.5))
    axis.vlines(
        [peak.mz for peak in peak_list],
        0,
        [100.0 * peak.intensity / maximum if maximum else 0.0 for peak in peak_list],
        color="#203864",
        linewidth=1.2,
    )
    axis.axvline(precursor_mz, color="#b22222", linestyle="--", linewidth=1, label=f"precursor {precursor_mz:.4f}")
    axis.set(
        title=title,
        xlabel="m/z",
        ylabel="Relative intensity (%)",
        ylim=(0, 108),
    )
    axis.grid(axis="y", alpha=0.2)
    axis.legend(loc="upper right")
    figure.tight_layout()
    figure.savefig(path, dpi=180)
    plt.close(figure)


def _draw_structure(smiles: str, path: Path) -> str:
    """Render a supplied SMILES to PNG and return its canonical SMILES."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
    except ImportError as error:
        raise RuntimeError("SMILES visualization requires RDKit; install requirements-structures.txt") from error
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    canonical_smiles = Chem.MolToSmiles(molecule)
    image = Draw.MolToImage(molecule, size=(900, 650))
    image.save(path)
    return canonical_smiles


def analyze_text_file(
    input_path: Path,
    output_dir: Path,
    top: int = 10,
    tolerance_ppm: float | None = None,
    smiles: str | None = None,
    allowed_elements: Iterable[str] | None = None,
    required_elements: Iterable[str] | None = None,
) -> list[Candidate]:
    """Write spectrum PNG, ranked-formula CSV, and a plain-text report.

    Explicit tolerance and element arguments override header metadata;
    otherwise use header values, then defaults (25 ppm and DEFAULT_ELEMENTS).
    An empty required-elements iterable clears header requirements. Supplied
    SMILES is rendered as a reference only and never used to rank candidates.
    Raise ValueError when precursor metadata or matching formulae are absent.
    """
    metadata, peaks = read_text_spectrum(input_path)
    if "precursor_mz" not in metadata:
        raise ValueError(f"{input_path}: add '# precursor_mz: ...' to the text header")
    precursor_mz = float(metadata["precursor_mz"])
    adduct = metadata.get("adduct", "proton").lower()
    tolerance_ppm = (
        float(metadata.get("tolerance_ppm", 25.0))
        if tolerance_ppm is None
        else tolerance_ppm
    )
    elements = parse_elements(
        allowed_elements
        if allowed_elements is not None
        else metadata.get("elements", ",".join(DEFAULT_ELEMENTS))
    )
    required_values = (
        required_elements if required_elements is not None
        else metadata.get("required_elements", "")
    )
    if not isinstance(required_values, str):
        required_values = tuple(required_values)
    required = parse_elements(required_values) if required_values else ()
    title = metadata.get("title", input_path.stem)
    candidates = rank_candidates(
        precursor_mz, peaks, adduct, tolerance_ppm, top, elements, required
    )
    if not candidates:
        raise ValueError(f"{input_path}: no formula matched the precursor within {tolerance_ppm:g} ppm")

    output_dir.mkdir(parents=True, exist_ok=True)
    _draw_spectrum(output_dir / "spectrum.png", title, precursor_mz, peaks)
    with (output_dir / "candidates.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow((
            "rank",
            "formula",
            "neutral_exact_mass",
            "ppm_error",
            "score",
            "fragment_mse",
            "matched_intensity_fraction",
            "explained_peaks",
        ))
        for index, candidate in enumerate(candidates, start=1):
            writer.writerow((
                index,
                candidate.formula,
                f"{candidate.exact_mass:.6f}",
                f"{candidate.ppm_error:.3f}",
                f"{candidate.score:.3f}",
                f"{candidate.fragment_mse:.6f}",
                f"{candidate.matched_intensity_fraction:.6f}",
                "; ".join(candidate.explained_peaks),
            ))

    report_lines = [
        f"Spectrum: {title}",
        f"Input: {input_path}",
        f"Precursor m/z: {precursor_mz:.5f} ({'[M]+' if adduct == 'none' else f'[M+{adduct}]+'})",
        f"Mass tolerance: {tolerance_ppm:g} ppm",
        f"Allowed elements: {', '.join(elements)}",
        f"Required elements: {', '.join(required) if required else 'none'}",
        "Formula suggestions (not unique structure identifications):",
    ]
    report_lines.extend(
        f"{index}. {candidate.formula}  exact mass {candidate.exact_mass:.5f}  "
        f"{candidate.ppm_error:+.2f} ppm  score {candidate.score:.3f}  "
        f"fragment MSE {candidate.fragment_mse:.3f}  "
        f"matched intensity {candidate.matched_intensity_fraction:.1%}"
        for index, candidate in enumerate(candidates, start=1)
    )
    if smiles:
        canonical_smiles = _draw_structure(smiles, output_dir / "structure.png")
        report_lines.extend((
            "Reference structure from the presentation (not inferred from MS alone):",
            f"Provided SMILES: {smiles}",
            f"Canonical SMILES: {canonical_smiles}",
        ))
    (output_dir / "report.txt").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print("\n".join(report_lines))
    print(f"Spectrum plot: {output_dir / 'spectrum.png'}")
    if smiles:
        print(f"Structure image: {output_dir / 'structure.png'}")
    print(f"Candidate table: {output_dir / 'candidates.csv'}")
    return candidates


def _parse_elements_argument(value: str) -> tuple[str, ...]:
    """Convert an argparse element list into validated element symbols."""
    try:
        return parse_elements(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error


DEMO_PRECURSOR_MZ = 257.16
DEMO_PEAKS = (
    Peak(113.10, 2),
    Peak(130.07, 2),
    Peak(145.08, 2),
    Peak(157.08, 2),
    Peak(169.08, 2),
    Peak(183.09, 8),
    Peak(185.11, 24),
    Peak(186.11, 6),
    Peak(197.11, 3),
    Peak(223.12, 4),
    Peak(225.14, 85),
    Peak(226.14, 10),
    Peak(228.13, 5),
    Peak(241.13, 8),
    Peak(253.17, 18),
    Peak(256.16, 38),
    Peak(257.16, 100),
    Peak(258.17, 6),
    Peak(271.18, 3),
    Peak(285.20, 20),
    Peak(286.20, 7),
    Peak(298.20, 2),
)


def main() -> int:
    """Parse CLI arguments and run demo, single-file, or batch analysis."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="text spectrum file (header metadata plus m/z,intensity rows)")
    parser.add_argument("--batch", type=Path, help="directory of .txt spectrum inputs to analyze")
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--smiles", help="optional known/reference SMILES to visualize for one input")
    parser.add_argument("--structures", type=Path, help="CSV mapping spectrum_file,smiles for batch reference structures")
    parser.add_argument("--precursor", type=float, help="precursor m/z")
    parser.add_argument("--adduct", choices=ADDUCT_MASSES, default="proton")
    parser.add_argument("--tolerance-ppm", type=float)
    parser.add_argument(
        "--elements",
        type=_parse_elements_argument,
        help="allowed precursor elements, comma-separated (default: C,H,N,O,S,P,Cl,Br)",
    )
    parser.add_argument(
        "--require-elements",
        type=_parse_elements_argument,
        default=(),
        help="elements that must occur at least once in each candidate formula",
    )
    parser.add_argument("--peak", action="append", type=_parse_peak, default=[], help="m/z:intensity")
    parser.add_argument("--peaks", type=Path, help="CSV file with mz,intensity columns")
    parser.add_argument("--top", type=int, default=10, help="number of formula candidates to show")
    parser.add_argument("--demo", action="store_true", help="run on the first spectrum in Data.pptx")
    args = parser.parse_args()
    cli_tolerance_ppm = (
        args.tolerance_ppm
        if args.tolerance_ppm is not None
        else (50.0 if args.demo else 25.0)
    )

    if args.input or args.batch:
        if args.demo or args.precursor is not None or args.peak or args.peaks:
            parser.error("file analysis cannot be combined with --demo, --precursor, or inline peaks")
        if args.input and args.batch:
            parser.error("choose either --input or --batch")
        if args.input and args.structures:
            parser.error("--structures is only used with --batch")
        if args.batch:
            if not args.batch.is_dir():
                parser.error(f"batch path is not a directory: {args.batch}")
            input_paths = sorted(args.batch.glob("*.txt"))
            if not input_paths:
                parser.error(f"no .txt spectrum files found in {args.batch}")
            if args.smiles:
                parser.error("use --structures, not --smiles, with --batch")
            structure_map: dict[str, str] = {}
            if args.structures:
                with args.structures.open(newline="", encoding="utf-8") as handle:
                    reader = csv.DictReader(handle)
                    if not reader.fieldnames or not {"spectrum_file", "smiles"} <= set(reader.fieldnames):
                        parser.error("structure mapping CSV must have spectrum_file,smiles columns")
                    structure_map = {
                        row["spectrum_file"]: row["smiles"]
                        for row in reader
                        if row.get("spectrum_file") and row.get("smiles")
                    }
            batch_summary: list[tuple[str, dict[str, str], list[Candidate], str | None]] = []
            for input_path in input_paths:
                metadata, _ = read_text_spectrum(input_path)
                candidates = analyze_text_file(
                    input_path,
                    args.output_dir / input_path.stem,
                    args.top,
                    args.tolerance_ppm,
                    structure_map.get(input_path.name),
                    args.elements,
                    args.require_elements or None,
                )
                batch_summary.append(
                    (input_path.name, metadata, candidates, structure_map.get(input_path.name))
                )
            args.output_dir.mkdir(parents=True, exist_ok=True)
            with (args.output_dir / "summary.csv").open(
                "w", newline="", encoding="utf-8"
            ) as handle:
                writer = csv.writer(handle)
                writer.writerow((
                    "spectrum_file",
                    "precursor_mz",
                    "adduct",
                    "top_formula",
                    "top_formula_ppm",
                    "top_score",
                    "top_fragment_mse",
                    "top_matched_intensity_fraction",
                    "other_top_formulas",
                    "allowed_elements",
                    "required_elements",
                    "reference_smiles",
                ))
                for filename, metadata, candidates, reference_smiles in batch_summary:
                    writer.writerow((
                        filename,
                        metadata.get("precursor_mz", ""),
                        metadata.get("adduct", "proton"),
                        candidates[0].formula,
                        f"{candidates[0].ppm_error:.3f}",
                        f"{candidates[0].score:.3f}",
                        f"{candidates[0].fragment_mse:.6f}",
                        f"{candidates[0].matched_intensity_fraction:.6f}",
                        "; ".join(candidate.formula for candidate in candidates[1:]),
                        ",".join(
                            parse_elements(
                                args.elements
                                if args.elements is not None
                                else metadata.get("elements", DEFAULT_ELEMENTS)
                            )
                        ),
                        ",".join(
                            parse_elements(
                                args.require_elements
                                if args.require_elements
                                else metadata.get("required_elements", "")
                            )
                        ) if args.require_elements or metadata.get("required_elements") else "",
                        reference_smiles or "",
                    ))
            print(f"Batch summary: {args.output_dir / 'summary.csv'}")
        else:
            analyze_text_file(
                args.input,
                args.output_dir,
                args.top,
                args.tolerance_ppm,
                args.smiles,
                args.elements,
                args.require_elements or None,
            )
        return 0

    if args.smiles or args.structures:
        parser.error("--smiles and --structures require --input or --batch")
    if args.demo:
        if args.precursor is not None or args.peaks or args.peak:
            parser.error("--demo cannot be combined with --precursor, --peaks, or --peak")
        precursor_mz, peaks = DEMO_PRECURSOR_MZ, DEMO_PEAKS
    else:
        if args.precursor is None:
            parser.error("--precursor is required unless --demo is used")
        if args.peaks and args.peak:
            parser.error("use either --peaks or --peak, not both")
        precursor_mz = args.precursor
        peaks = _read_csv(args.peaks) if args.peaks else args.peak

    candidates = rank_candidates(
        precursor_mz,
        peaks,
        args.adduct,
        cli_tolerance_ppm,
        args.top,
        args.elements if args.elements is not None else DEFAULT_ELEMENTS,
        args.require_elements,
    )
    if not candidates:
        print("No formula candidates matched the precursor within the requested tolerance.")
        return 1

    adduct_label = "[M]+" if args.adduct == "none" else f"[M+{args.adduct}]+"
    elements = args.elements if args.elements is not None else DEFAULT_ELEMENTS
    print(
        f"Precursor m/z {precursor_mz:.5f}, adduct {adduct_label}, "
        f"tolerance {cli_tolerance_ppm:g} ppm"
    )
    print(f"Allowed elements: {', '.join(elements)}")
    print("Rank  Formula          Neutral mass   Error (ppm)  Score    Frag MSE  Matched intensity  Explained peaks")
    for rank, candidate in enumerate(candidates, start=1):
        explanations = ", ".join(candidate.explained_peaks) or "none"
        print(
            f"{rank:>4}  {candidate.formula:<15} "
            f"{candidate.exact_mass:>12.5f}  {candidate.ppm_error:>11.2f}  "
            f"{candidate.score:>7.3f}  {candidate.fragment_mse:>8.3f}  "
            f"{candidate.matched_intensity_fraction:>16.1%}  {explanations}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
