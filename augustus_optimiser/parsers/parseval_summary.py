import re

from typing import NamedTuple
from typing import Optional
from typing import List, Mapping, Dict, Tuple, Iterator
from typing import Any, Union

from augustus_optimiser.higher import fmap


WHITESPACE = re.compile(r"\s+")
PERIODS = re.compile(r"\.\.+")
COLON_WHITESPACE = re.compile(r":\s+")


def stripper(handle: Iterator[str]) -> Iterator[str]:
    """ Iterator adapter that strips leading and trailing whitespace. """

    for line in handle:
        yield line.strip()

    return


def split_header_row(line: str) -> str:
    """ Split by colon then whitespace and return right hand side. """

    split = COLON_WHITESPACE.split(line, maxsplit=1)

    if len(split) != 2:
        raise ValueError(f"Could not split '{line}' into two columns.")

    return split[1]


def split_dotted_row(line: str) -> str:
    """ Split by 2 or more periods and return right hand side. """

    split = PERIODS.split(line, maxsplit=1)

    if not len(split) == 2:
        raise ValueError(f"Could not split '{line}' into two columns.")

    return split[1]


def split_dotted_row_perc(line):
    """ Split by multiple periods and return number excluding bracketed bits.

    e.g. lines of the form:
    Test..............100 (0.4%)

    should return 100
    """

    split = split_dotted_row(line)

    regex = r"(\d+) \(.*\)"
    m = re.match(regex, split)
    if m is not None:
        return m.groups()[0]
    else:
        return "0"


def split_dotted_row_bp(line):
    """ Split by multiple periods and return the number excluding 'bp' """

    split = split_dotted_row(line)

    regex = r"(\d+\.?\d*) bp"
    m = re.match(regex, split)
    if m is not None:
        return m.groups()[0]
    else:
        return "0"


def split_dotted_row_aa(line):
    """ Split by multiple periods and return the number excluding 'aa' """

    split = split_dotted_row(line)

    regex = r"(\d+\.?\d*) aa"
    m = re.match(regex, split)
    if m is not None:
        return m.groups()[0]
    else:
        return "0"


def split_table_line(line):
    """
    """

    split = COLON_WHITESPACE.split(line, maxsplit=1)

    if len(split) != 2:
        raise ValueError(f"Could not split {line} into two columns.")

    split_floats = WHITESPACE.split(split[1], maxsplit=2)
    if len(split_floats) != 3:
        raise ValueError(f"Could not split '{split[1]}' into three columns.")

    return tuple(split_floats)


def is_null(string: Any) -> bool:
    return (string == "--") or (string is None) or (string == "")


def error_if_null(string: Any) -> None:
    if is_null(string):
        raise ValueError(
            "Received NULL value (--) but don't support optional values here. "
            f"Line was {string}."
        )


def any_or_none(string: Any) -> Optional[Any]:
    if is_null(string):
        return None
    else:
        return string


def int_or_error(string: Any) -> int:
    error_if_null(string)
    return int(string)


def float_or_error(string: Any) -> float:
    error_if_null(string)
    return float(string)


def str_or_error(string: Any) -> str:
    error_if_null(string)
    return str(string)


def float_or_none(string: Any) -> Optional[float]:
    return fmap(float, any_or_none(string))


def int_or_none(string: Any) -> Optional[int]:
    return fmap(int, any_or_none(string))


def get_or_none(obj: Mapping[Any, Any], key: Any) -> Optional[Any]:
    return obj.get(key, None)


def parse_sequences(handle):
    out = list()

    for line in stripper(handle):
        if line == "":
            break
        else:
            out.append(line)

    return out


def parse_gene_loci(firstline: str, handle: Iterator[str]) -> Dict[str, str]:
    """ Parse gene loci section into dict of strings. """

    out = dict()

    assert firstline.startswith("Gene loci")
    out["total"] = split_dotted_row(firstline)

    for line in stripper(handle):
        if line == "":
            break
        elif line.startswith("shared"):
            assert "shared" not in out
            out["shared"] = split_dotted_row(line)
        elif line.startswith("unique to reference"):
            assert "ref_unique" not in out
            out["ref_unique"] = split_dotted_row(line)
        elif line.startswith("unique to prediction"):
            assert "pred_unique" not in out
            out["pred_unique"] = split_dotted_row(line)
        else:
            raise ValueError(
                f"Received unexpected line in parseval output. {line}"
            )
    return out


class GeneLoci(NamedTuple):

    total: int
    shared: int
    ref_unique: int
    pred_unique: int

    @classmethod
    def from_obj(cls, obj: Mapping[Any, Any]) -> "GeneLoci":
        """ Parse a dict-like object as the class. """

        total = int_or_error(get_or_none(obj, "total"))
        shared = int_or_error(get_or_none(obj, "shared"))
        ref_unique = int_or_error(get_or_none(obj, "ref_unique"))
        pred_unique = int_or_error(get_or_none(obj, "pred_unique"))
        return cls(total, shared, ref_unique, pred_unique)

    def to_obj(self) -> Dict[str, int]:
        return {
            "total": self.total,
            "shared": self.shared,
            "ref_unique": self.ref_unique,
            "pred_unique": self.pred_unique,
        }


def parse_annotations(handle: Iterator[str]) -> Dict[str, str]:
    """ Parse annotations section into dict of strings. """

    out: Dict[str, str] = dict()

    context = None

    for line in stripper(handle):
        if line == "":
            break

        elif line.startswith("genes"):
            assert "num_genes" not in out
            out["num_genes"] = split_dotted_row(line)
            context = "genes"

        elif line.startswith("transcripts"):
            assert "num_transcripts" not in out
            out["num_transcripts"] = split_dotted_row(line)
            context = "transcripts"

        elif line.startswith("average per locus"):
            if context == "genes":
                assert "avg_genes_per_locus" not in out
                out["avg_genes_per_locus"] = split_dotted_row(line)

            elif context == "transcripts":
                assert "avg_transcripts_per_locus" not in out
                out["avg_transcripts_per_locus"] = split_dotted_row(line)

            else:
                raise ValueError(
                    f"Cannot parse line {line} with context: {context}"
                )
        elif line.startswith("average per gene"):
            if context == "transcripts":
                assert "avg_transcripts_per_gene" not in out
                out["avg_transcripts_per_gene"] = split_dotted_row(line)

            else:
                raise ValueError(
                    f"Cannot parse line {line} with context: {context}"
                )

    return out


class Annotations(NamedTuple):

    num_genes: int
    num_transcripts: int
    avg_genes_per_locus: Optional[float]
    avg_transcripts_per_locus: Optional[float]
    avg_transcripts_per_gene: Optional[float]

    @classmethod
    def from_obj(cls, obj: Mapping[Any, Any]) -> "Annotations":
        """ Convert from dict-like object. """

        num_genes = int_or_error(get_or_none(obj, "num_genes"))

        num_transcripts = int_or_error(get_or_none(obj, "num_transcripts"))

        avg_genes_per_locus = float_or_none(
            get_or_none(obj, "avg_genes_per_locus")
        )

        avg_transcripts_per_locus = float_or_none(
            get_or_none(obj, "avg_transcripts_per_locus")
        )

        avg_transcripts_per_gene = float_or_none(
            get_or_none(obj, "avg_transcripts_per_gene")
        )

        return cls(
            num_genes,
            num_transcripts,
            avg_genes_per_locus,
            avg_transcripts_per_locus,
            avg_transcripts_per_gene
        )

    def to_obj(self) -> Dict[str, Union[int, Optional[float]]]:
        return {
            "num_genes": self.num_genes,
            "num_transcripts": self.num_transcripts,
            "avg_genes_per_locus": self.avg_genes_per_locus,
            "avg_transcripts_per_locus": self.avg_transcripts_per_locus,
            "avg_transcripts_per_gene": self.avg_transcripts_per_gene,
        }


def parse_total_comparisons_matches(
    firstline: str,
    handle: Iterator[str]
) -> Tuple[Dict[str, str], str]:
    """ Parse a total comparisons subsection. """

    out = dict()
    out["num_matches"] = split_dotted_row_perc(firstline)

    for line in stripper(handle):
        if (
            line == ""
            or line.startswith("perfect matches")
            or line.startswith("CDS structure matches")
            or line.startswith("exon structure matches")
            or line.startswith("UTR structure matches")
            or line.startswith("non-matches")
        ):
            break
        elif line.startswith("avg. length"):
            assert "avg_length_bp" not in out
            out["avg_length_bp"] = split_dotted_row_bp(line)
        elif line.startswith("avg. # refr exons"):
            assert "avg_num_ref_exons" not in out
            out["avg_num_ref_exons"] = split_dotted_row(line)
        elif line.startswith("avg. # pred exons"):
            assert "avg_num_pred_exons" not in out
            out["avg_num_pred_exons"] = split_dotted_row(line)
        elif line.startswith("avg. refr CDS length"):
            assert "avg_ref_cds_length_aa" not in out
            out["avg_ref_cds_length_aa"] = split_dotted_row_aa(line)
        elif line.startswith("avg. pred CDS length"):
            assert "avg_pred_cds_length_aa" not in out
            out["avg_pred_cds_length_aa"] = split_dotted_row_aa(line)
        else:
            print(out)
            raise ValueError(
                f"Received unexpected line in total comparisons block: {line}"
            )

    return out, line


class ComparisonMatches(NamedTuple):

    num_matches: int
    avg_length_bp: Optional[float]
    avg_num_ref_exons: Optional[float]
    avg_num_pred_exons: Optional[float]
    avg_ref_cds_length_aa: Optional[float]
    avg_pred_cds_length_aa: Optional[float]

    @classmethod
    def from_obj(cls, obj: Mapping[Any, Any]) -> "ComparisonMatches":
        """ Parse from dict-like object """

        num_matches = int_or_error(get_or_none(obj, "num_matches"))
        avg_length_bp = float_or_none(get_or_none(obj, "avg_length_bp"))
        avg_num_ref_exons = float_or_none(
            get_or_none(obj, "avg_num_ref_exons")
        )

        avg_num_pred_exons = float_or_none(
            get_or_none(obj, "avg_num_pred_exons")
        )

        avg_ref_cds_length_aa = float_or_none(
            get_or_none(obj, "avg_ref_cds_length_aa")
        )

        avg_pred_cds_length_aa = float_or_none(
            get_or_none(obj, "avg_pred_cds_length_aa")
        )
        return cls(
            num_matches,
            avg_length_bp,
            avg_num_ref_exons,
            avg_num_pred_exons,
            avg_ref_cds_length_aa,
            avg_pred_cds_length_aa,
        )

    def to_obj(self) -> Dict[str, Union[int, Optional[float]]]:
        """ """

        return {
            "num_matches": self.num_matches,
            "avg_length_bp": self.avg_length_bp,
            "avg_num_ref_exons": self.avg_num_ref_exons,
            "avg_num_pred_exons": self.avg_num_pred_exons,
            "avg_ref_cds_length_aa": self.avg_ref_cds_length_aa,
            "avg_pred_cds_length_aa": self.avg_pred_cds_length_aa,
        }


def choose_total_comparisons_block(
    line: str,
    handle: Iterator[str],
    out: Dict[str, Union[str, Dict[str, str]]],
) -> Optional[str]:
    """ Does the heavy lifting for parse_total_comparisons.

    This split seems to be necessary because it's not easy to switch between
    subblocks without context, and recursion is the easiest way to do
    that.
    """

    if line == "":
        return line

    elif line.startswith("perfect matches with mislabeled UTRs"):
        assert "num_perfect_matches_mislabeled_utr" not in out
        out["num_perfect_matches_mislabeled_utr"] = split_dotted_row_perc(line)
        last_line = None

    elif line.startswith("perfect matches"):
        assert "perfect_matches" not in out
        pm, last_line = parse_total_comparisons_matches(line, handle)
        out["perfect_matches"] = pm

    elif line.startswith("CDS structure matches"):
        assert "cds_matches" not in out
        cm, last_line = parse_total_comparisons_matches(line, handle)
        out["cds_matches"] = cm

    elif line.startswith("exon structure matches"):
        assert "exon_matches" not in out
        em, last_line = parse_total_comparisons_matches(line, handle)
        out["exon_matches"] = em

    elif line.startswith("UTR structure matches"):
        assert "utr_matches" not in out
        um, last_line = parse_total_comparisons_matches(line, handle)
        out["utr_matches"] = um

    elif line.startswith("non-matches"):
        assert "non_matches" not in out
        nm, last_line = parse_total_comparisons_matches(line, handle)
        out["non_matches"] = nm

    else:
        print(out)
        raise ValueError(
            f"Received unexpected line in total comparisons section: {line}"
        )

    # Could remove this recursion with a peeking iterator.
    if last_line is not None and last_line != "":
        last_line = choose_total_comparisons_block(last_line, handle, out)

    if last_line == "":
        return ""

    return None


def parse_total_comparisons(
    firstline: str,
    handle: Iterator[str],
) -> Dict[str, Union[str, Dict[str, str]]]:
    """ Parse the total comparisons block.

    Farms out the heavy lifting to `choose_total_comparisons_block`.
    """

    out: Dict[str, Union[str, Dict[str, str]]] = dict()

    assert firstline.startswith("Total comparisons")
    out["total_comparisons"] = split_dotted_row(firstline)

    for line in stripper(handle):
        if line == "":
            break

        last_line = choose_total_comparisons_block(line, handle, out)
        if last_line == "":
            break

    return out


class TotalComparisons(NamedTuple):

    total_comparisons: int
    perfect_matches: ComparisonMatches
    num_perfect_matches_mislabeled_utr: int
    cds_matches: ComparisonMatches
    exon_matches: ComparisonMatches
    utr_matches: ComparisonMatches
    non_matches: ComparisonMatches

    @classmethod
    def from_obj(cls, obj: Mapping[Any, Any]) -> "TotalComparisons":
        """ Parse dict like object """

        total_comparisons = int_or_error(get_or_none(obj, "total_comparisons"))
        num_perfect_matches_mislabeled_utr = int_or_error(
            get_or_none(obj, "num_perfect_matches_mislabeled_utr")
        )

        perfect_matches = ComparisonMatches.from_obj(
            obj.get("perfect_matches", {})
        )

        cds_matches = ComparisonMatches.from_obj(
            obj.get("cds_matches", {})
        )

        exon_matches = ComparisonMatches.from_obj(
            obj.get("exon_matches", {})
        )

        utr_matches = ComparisonMatches.from_obj(
            obj.get("utr_matches", {})
        )

        non_matches = ComparisonMatches.from_obj(
            obj.get("non_matches", {})
        )

        return cls(
            total_comparisons,
            perfect_matches,
            num_perfect_matches_mislabeled_utr,
            cds_matches,
            exon_matches,
            utr_matches,
            non_matches,
        )

    def to_obj(
        self
    ) -> Dict[str, Union[int, Dict[str, Union[int, Optional[float]]]]]:
        return {
            "total_comparisons": self.total_comparisons,
            "perfect_matches": self.perfect_matches.to_obj(),
            "num_perfect_matches_mislabeled_utr": self.num_perfect_matches_mislabeled_utr,  # noqa
            "cds_matches": self.cds_matches.to_obj(),
            "exon_matches": self.exon_matches.to_obj(),
            "utr_matches": self.utr_matches.to_obj(),
            "non_matches": self.non_matches.to_obj(),
        }


def parse_structure(handle: Iterator[str]) -> Dict[str, str]:  # noqa
    """ Parse the CDS/exon/utr structure sections. """

    out: Dict[str, str] = dict()

    for line in stripper(handle):
        if line == "":
            break

        elif line.startswith("reference"):
            assert "num_ref_segments" not in out
            out["num_ref_segments"] = split_dotted_row(line)

        elif line.startswith("prediction"):
            assert "num_pred_segments" not in out
            out["num_pred_segments"] = split_dotted_row(line)

        elif line.startswith("match prediction"):
            assert "num_ref_match_pred" not in out
            out["num_ref_match_pred"] = split_dotted_row_perc(line)

        elif line.startswith("don't match prediction"):
            assert "num_ref_dont_match_pred" not in out
            out["num_ref_dont_match_pred"] = split_dotted_row_perc(line)

        elif line.startswith("match reference"):
            assert "num_pred_match_ref" not in out
            out["num_pred_match_ref"] = split_dotted_row_perc(line)

        elif line.startswith("don't match reference"):
            assert "num_pred_dont_match_ref" not in out
            out["num_pred_dont_match_ref"] = split_dotted_row_perc(line)

        elif line.startswith("Sensitivity"):
            assert "sensitivity" not in out
            out["sensitivity"] = split_dotted_row(line)

        elif line.startswith("Specificity"):
            assert "specificity" not in out
            out["specificity"] = split_dotted_row(line)

        elif line.startswith("F1 Score"):
            assert "f1_score" not in out
            out["f1_score"] = split_dotted_row(line)

        elif line.startswith("Annotation edit distance"):
            assert "edit_distance" not in out
            out["edit_distance"] = split_dotted_row(line)

        else:
            raise ValueError(f"Unexpected line in structure block: {line}")

    return out


class StructureComparisons(NamedTuple):

    num_ref_segments: int
    num_pred_segments: int
    num_ref_match_pred: int
    num_ref_dont_match_pred: int
    num_pred_match_ref: int
    num_pred_dont_match_ref: int
    sensitivity: Optional[float]
    specificity: Optional[float]
    f1_score: Optional[float]
    edit_distance: Optional[float]

    @classmethod
    def from_obj(cls, obj: Mapping[Any, Any]) -> "StructureComparisons":
        """ Parse from a dict-like object. """

        num_ref_segments = int_or_error(get_or_none(obj, "num_ref_segments"))
        num_pred_segments = int_or_error(get_or_none(obj, "num_pred_segments"))

        num_ref_match_pred = int_or_error(
            get_or_none(obj, "num_ref_match_pred")
        )

        num_ref_dont_match_pred = int_or_error(
            get_or_none(obj, "num_ref_dont_match_pred")
        )

        num_pred_match_ref = int_or_error(
            get_or_none(obj, "num_pred_match_ref")
        )

        num_pred_dont_match_ref = int_or_error(
            get_or_none(obj, "num_pred_dont_match_ref")
        )

        sensitivity = float_or_none(get_or_none(obj, "sensitivity"))
        specificity = float_or_none(get_or_none(obj, "specificity"))
        f1_score = float_or_none(get_or_none(obj, "f1_score"))
        edit_distance = float_or_none(get_or_none(obj, "edit_distance"))

        return cls(
            num_ref_segments,
            num_pred_segments,
            num_ref_match_pred,
            num_ref_dont_match_pred,
            num_pred_match_ref,
            num_pred_dont_match_ref,
            sensitivity,
            specificity,
            f1_score,
            edit_distance,
        )

    def to_obj(self) -> Dict["str", Union[int, Optional[float]]]:
        return {
            "num_ref_segments": self.num_ref_segments,
            "num_pred_segments": self.num_pred_segments,
            "num_ref_match_pred": self.num_ref_match_pred,
            "num_ref_dont_match_pred": self.num_ref_dont_match_pred,
            "num_pred_match_ref": self.num_pred_match_ref,
            "num_pred_dont_match_ref": self.num_pred_dont_match_ref,
            "sensitivity": self.sensitivity,
            "specificity": self.specificity,
            "f1_score": self.f1_score,
            "edit_distance": self.edit_distance,
        }


def parse_nucleotide_table(handle):
    """ Parse the final nucleotide-level statistics table.
    This comes in three columns, so has to be parsed together.
    """

    cds = dict()
    utr = dict()
    overall = dict()

    def add_match(key, c, u, o, cds, utr, overall):
        assert key not in cds
        assert key not in utr
        assert key not in overall

        cds[key] = c
        utr[key] = u
        overall[key] = o
        return

    for line in stripper(handle):
        if line == "":
            break
        elif line.startswith("Matching coefficient"):
            c, u, o = split_table_line(line)
            add_match("matching_coefficient", c, u, o, cds, utr, overall)
        elif line.startswith("Correlation coefficient"):
            c, u, o = split_table_line(line)
            add_match("correlation_coefficient", c, u, o, cds, utr, overall)
        elif line.startswith("Sensitivity"):
            c, u, o = split_table_line(line)
            add_match("sensitivity", c, u, o, cds, utr, overall)
        elif line.startswith("Specificity"):
            c, u, o = split_table_line(line)
            add_match("specificity", c, u, o, cds, utr, overall)
        elif line.startswith("F1 Score"):
            c, u, o = split_table_line(line)
            add_match("f1_score", c, u, o, cds, utr, overall)
        elif line.startswith("Annotation edit distance"):
            c, u, o = split_table_line(line)
            add_match("edit_distance", c, u, o, cds, utr, overall)
        else:
            raise ValueError(
                f"Encountered unexpected line in nucleotide block: {line}"
            )

    return cds, utr, overall


class NucleotideComparison(NamedTuple):

    matching_coefficient: Optional[float]
    correlation_coefficient: Optional[float]
    sensitivity: Optional[float]
    specificity: Optional[float]
    f1_score: Optional[float]
    edit_distance: Optional[float]

    @classmethod
    def from_obj(cls, obj: Mapping[Any, Any]) -> "NucleotideComparison":
        matching_coefficient = float_or_none(
            get_or_none(obj, "matching_coefficient")
        )

        correlation_coefficient = float_or_none(
            get_or_none(obj, "correlation_coefficient")
        )

        sensitivity = float_or_none(get_or_none(obj, "sensitivity"))
        specificity = float_or_none(get_or_none(obj, "specificity"))
        f1_score = float_or_none(get_or_none(obj, "f1_score"))
        edit_distance = float_or_none(get_or_none(obj, "edit_distance"))
        return cls(
            matching_coefficient,
            correlation_coefficient,
            sensitivity,
            specificity,
            f1_score,
            edit_distance,
        )

    def to_obj(self) -> Dict[str, Optional[float]]:
        return {
            "matching_coefficient": self.matching_coefficient,
            "correlation_coefficient": self.correlation_coefficient,
            "sensitivity": self.sensitivity,
            "specificity": self.specificity,
            "f1_score": self.f1_score,
            "edit_distance": self.edit_distance,
        }


def parse_parseval_summary(handle: Iterator[str]) -> Dict[str, Any]:  # noqa
    """ Parse the parseval summary file.

    Parseval needs to be run in summary text mode '-s --outformat text' to
    prevent the per-gene gff3 comparisons.
    """

    out: Dict[str, Any] = dict()

    for line in stripper(handle):

        if line.startswith("="):
            continue

        elif line == "":
            continue

        elif line.startswith("Started:"):
            assert "started" not in out
            out["started"] = split_header_row(line)

        # Note that the colon is important.
        elif line.startswith("Reference annotations:"):
            assert "reference" not in out
            out["reference"] = split_header_row(line)

        # Note that the colon is important.
        elif line.startswith("Prediction annotations:"):
            assert "prediction" not in out
            out["prediction"] = split_header_row(line)

        elif line.startswith("Executing command"):
            assert "cmd" not in out
            out["cmd"] = split_header_row(line)

        elif line.startswith("Sequences compared"):
            assert "sequences" not in out
            out["sequences"] = parse_sequences(handle)

        elif line.startswith("Gene loci"):
            assert "gene_loci" not in out
            out["gene_loci"] = parse_gene_loci(line, handle)

        elif line.startswith("Reference annotations"):
            assert "ref_annotations" not in out
            out["ref_annotations"] = parse_annotations(handle)

        elif line.startswith("Prediction annotations"):
            assert "pred_annotations" not in out
            out["pred_annotations"] = parse_annotations(handle)

        elif line.startswith("Total comparisons"):
            assert "total_comparisons" not in out
            out["total_comparisons"] = parse_total_comparisons(line, handle)

        elif line.startswith("CDS structure comparison"):
            assert "cds_structure" not in out
            out["cds_structure"] = parse_structure(handle)

        elif line.startswith("Exon structure comparison"):
            assert "exon_structure" not in out
            out["exon_structure"] = parse_structure(handle)

        elif line.startswith("UTR structure comparison"):
            assert "utr_structure" not in out
            out["utr_structure"] = parse_structure(handle)

        elif line.startswith("Nucleotide-level comparison"):
            cds, utr, overall = parse_nucleotide_table(handle)

            assert "cds_nucleotide_comparisons" not in out
            out["cds_nucleotide_comparisons"] = cds

            assert "utr_nucleotide_comparisons" not in out
            out["utr_nucleotide_comparisons"] = utr

            assert "overall_nucleotide_comparisons" not in out
            out["overall_nucleotide_comparisons"] = overall

        else:
            print(out)
            raise ValueError(f"unexpected {line}")

    return out


class ParsEvalSummary(NamedTuple):

    started: str
    reference: str
    prediction: str
    cmd: str
    sequences: List[str]
    gene_loci: GeneLoci
    ref_annotations: Annotations
    pred_annotations: Annotations
    total_comparisons: TotalComparisons
    cds_structure: StructureComparisons
    exon_structure: StructureComparisons
    utr_structure: StructureComparisons
    cds_nucleotide_comparison: NucleotideComparison
    utr_nucleotide_comparison: NucleotideComparison
    overall_nucleotide_comparison: NucleotideComparison

    @classmethod
    def from_obj(cls, obj: Mapping[Any, Any]) -> "ParsEvalSummary":
        """ """

        started = str_or_error(get_or_none(obj, "started"))
        reference = str_or_error(get_or_none(obj, "reference"))
        prediction = str_or_error(get_or_none(obj, "prediction"))
        cmd = str_or_error(get_or_none(obj, "cmd"))
        sequences = [str(o) for o in obj.get("sequences", [])]
        gene_loci = GeneLoci.from_obj(obj.get("gene_loci", {}))
        ref_annotations = Annotations.from_obj(obj.get("ref_annotations", {}))

        pred_annotations = Annotations.from_obj(
            obj.get("pred_annotations", {})
        )

        total_comparisons = TotalComparisons.from_obj(
            obj.get("total_comparisons", {})
        )

        cds_structure = StructureComparisons.from_obj(
            obj.get("cds_structure", {})
        )

        exon_structure = StructureComparisons.from_obj(
            obj.get("exon_structure", {})
        )

        utr_structure = StructureComparisons.from_obj(
            obj.get("utr_structure", {})
        )

        cds_nucleotide_comparison = NucleotideComparison.from_obj(
            obj.get("cds_nucleotide_comparison", {})
        )

        utr_nucleotide_comparison = NucleotideComparison.from_obj(
            obj.get("utr_nucleotide_comparison", {})
        )

        overall_nucleotide_comparison = NucleotideComparison.from_obj(
            obj.get("overall_nucleotide_comparison", {})
        )

        return cls(
            started,
            reference,
            prediction,
            cmd,
            sequences,
            gene_loci,
            ref_annotations,
            pred_annotations,
            total_comparisons,
            cds_structure,
            exon_structure,
            utr_structure,
            cds_nucleotide_comparison,
            utr_nucleotide_comparison,
            overall_nucleotide_comparison,
        )

    @classmethod
    def parse(cls, handle: Iterator[str]) -> "ParsEvalSummary":
        obj = parse_parseval_summary(handle)
        try:
            self = cls.from_obj(obj)
        except Exception as e:
            print(obj)
            raise e

        return self
