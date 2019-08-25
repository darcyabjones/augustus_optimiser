import re

from typing import NamedTuple
from typing import Optional
from typing import List, Sequence, Tuple

WHITESPACE = re.compile(r"\s+")
PERIODS = re.compile(r"\.\.+")
COLON_WHITESPACE = re.compile(r":\s+")


def skip_blank_lines(handle: Sequence[str]):
    for line in handle:
        stripped = line.strip()
        if stripped != "" and not stripped.startswith("="):
            yield stripped
    return


def float_or_none(string: str) -> Optional[float]:
    try:
        return float(string)
    except ValueError:
        return None


def match_next_line(handle: Sequence[str], expected: str) -> None:
    stripped = next(handle).strip()
    match = re.fullmatch(expected, stripped)
    if match is None:
        raise ValueError(
            f"The line {stripped} didn't match what we expected {expected}."
        )
    return


def parse_split_line(handle: Sequence[str], expected: str) -> str:

    stripped = next(handle).strip()
    split = PERIODS.split(stripped, maxsplit=1)

    if not len(split) == 2:
        raise ValueError("Could not split line into two.")

    match = re.fullmatch(expected, split[0])
    if match is None:
        raise ValueError(f"Expected key to be {expected} got {split[0]}")

    return split[1]


def parse_int_line(string: str, expected_key: str) -> int:
    val = parse_split_line(string, expected_key)
    split_val = val.split()
    return int(split_val[0])


def parse_float_line(string: str, expected_key: str) -> float:
    val = parse_split_line(string, expected_key)
    split_val = val.split()
    return float(split_val[0])


class GeneLoci(NamedTuple):

    total: int
    shared: int
    unique_ref: int
    unique_pred: int

    @classmethod
    def parse(cls, handle: Sequence[str]) -> "GeneLoci":
        """

        Examples
        >>> string = '''
        ...  Gene loci................................14
        ...    shared.................................11
        ...    unique to reference....................1
        ...    unique to prediction...................2
        ... '''
        >>> GeneLoci.from_str(skip_blank_lines(string.split('\n')))
        GeneLoci(total=14, shared=11, unique_ref=1, unique_pred=2)
        """

        return cls(
            parse_int_line(handle, "Gene loci"),
            parse_int_line(handle, "shared"),
            parse_int_line(handle, "unique to reference"),
            parse_int_line(handle, "unique to prediction")
        )


class Comparisons(NamedTuple):

    ncomparisons: int
    perfect_matches: int
    perfect_matches_mislabeled_utrs: int
    cds_matches: int
    cds_avg_length_bp: float
    cds_avg_num_ref_exons: float
    cds_avg_num_pred_exons: float
    cds_avg_length_ref_cds_aa: float
    cds_avg_length_pred_cds_aa: float
    exon_matches: int
    utr_matches: int
    non_matches: int
    non_matches_avg_length_bp: float
    non_matches_avg_num_ref_exons: float
    non_matches_avg_num_pred_exons: float
    non_matches_avg_length_ref_cds_aa: float
    non_matches_avg_length_pred_cds_aa: float

    @classmethod
    def parse(cls, handle: Sequence[str]) -> "Comparisons":
        """
        >>> string = '''
        ...  Total comparisons........................11
        ...    perfect matches........................0 (0.0%)
        ...    perfect matches with mislabeled UTRs...0 (0.0%)
        ...    CDS structure matches..................4 (36.4%)
        ...      avg. length..........................3781.50 bp
        ...      avg. # refr exons....................7.25
        ...      avg. # pred exons....................7.25
        ...      avg. refr CDS length.................275.00 aa
        ...      avg. pred CDS length.................275.00 aa
        ...    exon structure matches.................0 (0.0%)
        ...    UTR structure matches..................0 (0.0%)
        ...    non-matches............................7 (63.6%)
        ...      avg. length..........................2807.14 bp
        ...      avg. # refr exons....................4.14
        ...      avg. # pred exons....................2.43
        ...      avg. refr CDS length.................281.86 aa
        ...      avg. pred CDS length.................384.86 aa
        ... '''
        >>> Comparisons.parse(skip_blank_lines(string.split('\n')))
        Comparisons(ncomparisons=11,
                    perfect_matches=0,
                    perfect_matches_mislabeled_utrs=0,
                    cds_matches=4,
                    cds_avg_length_bp=3781.5,
                    cds_avg_num_ref_exons=7.25,
                    cds_avg_num_pred_exons=7.25,
                    cds_avg_length_ref_cds_aa=275.0,
                    cds_avg_length_pred_cds_aa=275.0,
                    exon_matches=0,
                    utr_matches=0,
                    non_matches=7,
                    non_matches_avg_length_bp=2807.14,
                    non_matches_avg_num_ref_exons=4.14,
                    non_matches_avg_num_pred_exons=2.43,
                    non_matches_avg_length_ref_cds_aa=281.86,
                    non_matches_avg_length_pred_cds_aa=384.86)
        """

        return cls(
            parse_int_line(handle, "Total comparisons"),
            parse_int_line(handle, "perfect matches"),
            parse_int_line(handle, "perfect matches with mislabeled UTRs"),
            parse_int_line(handle, "CDS structure matches"),
            parse_float_line(handle, "avg. length"),
            parse_float_line(handle, "avg. # refr exons"),
            parse_float_line(handle, "avg. # pred exons"),
            parse_float_line(handle, "avg. refr CDS length"),
            parse_float_line(handle, "avg. pred CDS length"),
            parse_int_line(handle, "exon structure matches"),
            parse_int_line(handle, "UTR structure matches"),
            parse_int_line(handle, "non-matches"),
            parse_float_line(handle, "avg. length"),
            parse_float_line(handle, "avg. # refr exons"),
            parse_float_line(handle, "avg. # pred exons"),
            parse_float_line(handle, "avg. refr CDS length"),
            parse_float_line(handle, "avg. pred CDS length"),
        )


class StructureComparisons(NamedTuple):

    ref_segments: int
    ref_match_prediction: int
    ref_not_match_prediction: int
    pred_segments: int
    pred_match_reference: int
    pred_not_match_reference: int
    sensitivity: float
    specificity: float
    f1_score: float
    edit_distance: float

    @classmethod
    def parse(cls, handle: Sequence[str]) -> "StructureComparisons":
        """

        Examples:
        >>> string = '''
        ... reference CDS segments.................56
        ... match prediction.....................34 (60.7%)
        ...   don't match prediction...............22 (39.3%)
        ... prediction CDS segments................44
        ...   match reference......................34 (77.3%)
        ...   don't match reference................10 (22.7%)
        ... Sensitivity............................0.607
        ... Specificity............................0.773
        ... F1 Score...............................0.680
        ... Annotation edit distance...............0.310
        ... '''
        >>> StructureComparisons.parse(skip_blank_lines(string.split('\n')))
        StructureComparisons(ref_segments=56,
                             ref_match_prediction=34,
                             ref_not_match_prediction=22,
                             pred_segments=44,
                             pred_match_reference=34,
                             pred_not_match_reference=10,
                             sensitivity=0.607,
                             specificity=0.773,
                             f1_score=0.680,
                             edit_distance=0.310)
        """
        return cls(
            parse_int_line(handle, "reference .+"),
            parse_int_line(handle, "match prediction"),
            parse_int_line(handle, "don't match prediction"),
            parse_int_line(handle, "prediction .+"),
            parse_int_line(handle, "match reference"),
            parse_int_line(handle, "don't match reference"),
            parse_float_line(handle, "Sensitivity"),
            parse_float_line(handle, "Specificity"),
            parse_float_line(handle, "F1 Score"),
            parse_float_line(handle, "Annotation edit distance"),
        )


def split_table_line(
    handle: Sequence[str],
    expected: str
) -> Tuple[float, float, Optional[float]]:
    """
    """

    stripped = next(handle).strip()
    split = COLON_WHITESPACE.split(stripped, maxsplit=1)

    if len(split) != 2:
        raise ValueError(f"Could not split {stripped} into two columns.")

    match = re.fullmatch(expected, split[0])
    if match is None:
        raise ValueError(f"Expected key to be {expected} got {split[0]}")

    split_floats = WHITESPACE.split(split[1], maxsplit=2)
    if len(split_floats) != 3:
        raise ValueError(f"Could not split '{split[1]}' into three columns.")

    cds = float(split_floats[0])
    utrs = float(split_floats[1])
    overall = float_or_none(split_floats[2])
    return cds, utrs, overall


def split_metapar_line(handle: Sequence[str], expected: str) -> str:
    """
    """

    stripped = next(handle).strip()
    split = COLON_WHITESPACE.split(stripped, maxsplit=1)

    if len(split) != 2:
        raise ValueError(f"Could not split ''{stripped}'' into two columns.")

    match = re.fullmatch(expected, split[0])
    if match is None:
        raise ValueError(f"Expected key to be {expected} got {split[0]}")

    return split[1]


class NucleotideComparison(NamedTuple):

    matching_coefficient: float
    correlation_coefficient: float
    sensitivity: float
    specificity: float
    f1_score: float
    edit_distance: float

    @classmethod
    def parse(
        cls,
        handle: Sequence[str]
    ) -> Tuple["NucleotideComparison", "NucleotideComparison", float]:
        """
        >>> string = '''
        ...  Nucleotide-level comparison      CDS          UTRs         Overall
        ...    Matching coefficient:          0.904        0.902        0.840
        ...    Correlation coefficient:       0.780        0.536        --
        ...    Sensitivity:                   0.937        0.523        --
        ...    Specificity:                   0.759        0.666        --
        ...    F1 Score:                      0.838        0.586        --
        ...    Annotation edit distance:      0.152        0.405        --
        ... '''
        >>> c, u, o = NucleotideComparison.parse(
        ...     skip_blank_lines(string.split('\n'))
        ... )
        >>> c
        NucleotideComparison(matching_coefficient=0.904,
                             correlation_coefficient=0.78,
                             sensitivity=0.937,
                             specificity=0.759,
                             f1_score=0.838,
                             edit_distance=0.152)
        >>> o
        0.840
        """

        match_next_line(
            handle,
            r"Nucleotide-level comparison\s+CDS\s+UTRs\s+Overall"
        )
        mc_cds, mc_utr, mc_overall = split_table_line(
            handle,
            "Matching coefficient"
        )
        cc_cds, cc_utr, cc_overall = split_table_line(
            handle,
            "Correlation coefficient"
        )
        se_cds, se_utr, se_overall = split_table_line(handle, "Sensitivity")
        sp_cds, sp_utr, sp_overall = split_table_line(handle, "Specificity")
        f1_cds, f1_utr, f1_overall = split_table_line(handle, "F1 Score")
        ae_cds, ae_utr, ae_overall = split_table_line(
            handle,
            "Annotation edit distance"
        )

        cds = cls(mc_cds, cc_cds, se_cds, sp_cds, f1_cds, ae_cds)
        utr = cls(mc_utr, cc_utr, se_utr, sp_utr, f1_utr, ae_utr)

        if mc_overall is None:
            raise ValueError("Expected a Matching coefficient for Overall.")
        return cds, utr, mc_overall


class Annotations(NamedTuple):

    ngenes: int
    average_ngenes_per_locus: float
    ntranscripts: int
    average_ntranscripts_per_locus: float
    average_ntranscripts_per_gene: float

    @classmethod
    def parse(cls, handle: Sequence[str]) -> "Annotations":
        return cls(
            parse_int_line(handle, "genes"),
            parse_float_line(handle, "average per locus"),
            parse_int_line(handle, "transcripts"),
            parse_float_line(handle, "average per locus"),
            parse_float_line(handle, "average per gene")
        )


"""
    genes..................................12
      average per locus....................0.857
    transcripts............................12
      average per locus....................0.857
      average per gene.....................1.000

Annotations.parse(skip_blank_lines(string.split('\n')))
"""


def collect_until_blankline(handle: Sequence[str]) -> List[str]:
    out = []
    for line in handle:
        stripped = line.strip()
        if stripped == "":
            break
        else:
            out.append(stripped)
    return out


def parse_metapars(handle: Sequence[str]) -> Tuple[str, str, str]:
    st = split_metapar_line(handle, "Started")
    ref = split_metapar_line(handle, "Reference annotations")
    pred = split_metapar_line(handle, "Prediction annotations")
    _ = split_metapar_line(handle, "Executing command")
    return st, ref, pred


def parse_sequences(handle: Sequence[str]) -> List[str]:
    match_next_line(handle, "Sequences compared")
    return collect_until_blankline(handle)


class ParseEvalSummary(NamedTuple):

    started: str
    reference: str
    prediction: str
    sequences: List[str]
    gene_loci: GeneLoci
    ref_annotations: Annotations
    pred_annotations: Annotations
    total_comparisons: Comparisons
    cds_structure_comparison: StructureComparisons
    exon_structure_comparison: StructureComparisons
    utr_structure_comparison: StructureComparisons
    cds_nucleotide_comparison: NucleotideComparison
    utr_nucleotide_comparison: NucleotideComparison
    overall_nucleotide_matching_coefficient: float

    @classmethod
    def parse(cls, handle: Sequence[str]) -> "ParseEvalSummary":
        """
        """

        iterator = skip_blank_lines(handle)
        started, reference, prediction = parse_metapars(iterator)

        # Being cheeky here.
        # This is the only part that needs blank lines to delimit it.
        line = next(handle)
        assert line.strip() == ""
        sequences = parse_sequences(handle)

        gene_loci = GeneLoci.parse(iterator)

        _ = match_next_line(iterator, "Reference annotations")
        ref_annotations = Annotations.parse(iterator)

        _ = match_next_line(iterator, "Prediction annotations")
        pred_annotations = Annotations.parse(iterator)

        total_comparisons = Comparisons.parse(iterator)

        _ = match_next_line(iterator, "CDS structure comparison")
        cds_structure = StructureComparisons.parse(iterator)

        _ = match_next_line(iterator, "Exon structure comparison")
        exon_structure = StructureComparisons.parse(iterator)

        _ = match_next_line(iterator, "UTR structure comparison")
        utr_structure = StructureComparisons.parse(iterator)

        (cds_nuc_comp,
         utr_nuc_comp,
         overall_nucl) = NucleotideComparison.parse(iterator)

        return cls(
            started,
            reference,
            prediction,
            sequences,
            gene_loci,
            ref_annotations,
            pred_annotations,
            total_comparisons,
            cds_structure,
            exon_structure,
            utr_structure,
            cds_nuc_comp,
            utr_nuc_comp,
            overall_nucl
        )


class AnnotationStatistics(NamedTuple):

    num_gene_complete_match: int

    num_cds_str_match: int
    num_exon_str_match: int
    num_utr_str_match: int

    cds_num_true_pos: int
    cds_num_false_pos: int
    cds_num_false_neg: int
    cds_sens: float
    cds_spec: float
    cds_f1: float
    cds_edit: float

    exon_num_true_pos: int
    exon_num_false_pos: int
    exon_num_false_neg: int
    exon_sens: float
    exon_spec: float
    exon_f1: float
    exon_edit: float

    utr_num_true_pos: int
    utr_num_false_pos: int
    utr_num_false_neg: int
    utr_sens: float
    utr_spec: float
    utr_f1: float
    utr_edit: float

    cds_nucl_match_coef: float
    cds_nucl_corr_coef: float
    cds_nucl_sens: float
    cds_nucl_spec: float
    cds_nucl_f1: float
    cds_nucl_edit: float

    utr_nucl_match_coef: float
    utr_nucl_corr_coef: float
    utr_nucl_sens: float
    utr_nucl_spec: float
    utr_nucl_f1: float
    utr_nucl_edit: float

    overall_match_coef: float

    @staticmethod
    def columns() -> List[str]:
        return [
            "num_gene_complete_match",
            "num_cds_str_match",
            "num_exon_str_match",
            "num_utr_str_match",
            "cds_num_true_pos",
            "cds_num_false_pos",
            "cds_num_false_neg",
            "cds_sens",
            "cds_spec",
            "cds_f1",
            "cds_edit",
            "exon_num_true_pos",
            "exon_num_false_pos",
            "exon_num_false_neg",
            "exon_sens",
            "exon_spec",
            "exon_f1",
            "exon_edit",
            "utr_num_true_pos",
            "utr_num_false_pos",
            "utr_num_false_neg",
            "utr_sens",
            "utr_spec",
            "utr_f1",
            "utr_edit",
            "cds_nucl_match_coef",
            "cds_nucl_corr_coef",
            "cds_nucl_sens",
            "cds_nucl_spec",
            "cds_nucl_f1",
            "cds_nucl_edit",
            "utr_nucl_match_coef",
            "utr_nucl_corr_coef",
            "utr_nucl_sens",
            "utr_nucl_spec",
            "utr_nucl_f1",
            "utr_nucl_edit",
            "overall_match_coef",
        ]

    @classmethod
    def header(cls) -> str:
        return "\t".join(map(str, cls.columns()))

    def __str__(self) -> str:
        row = [
            str(getattr(self, c))
            for c
            in self.columns()
        ]

        return "\t".join(row)

    @classmethod
    def from_parse_eval_summary(
        cls,
        summary: ParseEvalSummary
    ) -> "AnnotationStatistics":
        return cls(
            summary.total_comparisons.perfect_matches,
            summary.total_comparisons.cds_matches,
            summary.total_comparisons.exon_matches,
            summary.total_comparisons.utr_matches,
            summary.cds_structure_comparison.pred_match_reference,
            summary.cds_structure_comparison.pred_not_match_reference,
            summary.cds_structure_comparison.ref_not_match_prediction,
            summary.cds_structure_comparison.sensitivity,
            summary.cds_structure_comparison.specificity,
            summary.cds_structure_comparison.f1_score,
            summary.cds_structure_comparison.edit_distance,
            summary.exon_structure_comparison.pred_match_reference,
            summary.exon_structure_comparison.pred_not_match_reference,
            summary.exon_structure_comparison.ref_not_match_prediction,
            summary.exon_structure_comparison.sensitivity,
            summary.exon_structure_comparison.specificity,
            summary.exon_structure_comparison.f1_score,
            summary.exon_structure_comparison.edit_distance,
            summary.utr_structure_comparison.pred_match_reference,
            summary.utr_structure_comparison.pred_not_match_reference,
            summary.utr_structure_comparison.ref_not_match_prediction,
            summary.utr_structure_comparison.sensitivity,
            summary.utr_structure_comparison.specificity,
            summary.utr_structure_comparison.f1_score,
            summary.utr_structure_comparison.edit_distance,
            summary.cds_nucleotide_comparison.matching_coefficient,
            summary.cds_nucleotide_comparison.correlation_coefficient,
            summary.cds_nucleotide_comparison.sensitivity,
            summary.cds_nucleotide_comparison.specificity,
            summary.cds_nucleotide_comparison.f1_score,
            summary.cds_nucleotide_comparison.edit_distance,
            summary.utr_nucleotide_comparison.matching_coefficient,
            summary.utr_nucleotide_comparison.correlation_coefficient,
            summary.utr_nucleotide_comparison.sensitivity,
            summary.utr_nucleotide_comparison.specificity,
            summary.utr_nucleotide_comparison.f1_score,
            summary.utr_nucleotide_comparison.edit_distance,
            summary.overall_nucleotide_matching_coefficient
        )
