from typing import NamedTuple
from typing import Any, Union, Optional
from typing import List, Mapping, Dict

from ao.higher import fmap
from ao.hints import HintConfig
from ao.parsers.parseval_summary import ParsEvalSummary


class AnnotationStatistics(NamedTuple):

    num_gene_complete_match: int

    num_cds_str_match: int
    num_exon_str_match: int
    num_utr_str_match: int

    cds_num_true_pos: int
    cds_num_false_pos: int
    cds_num_false_neg: int
    cds_sens: Optional[float]
    cds_spec: Optional[float]
    cds_f1: Optional[float]
    cds_edit: Optional[float]

    exon_num_true_pos: int
    exon_num_false_pos: int
    exon_num_false_neg: int
    exon_sens: Optional[float]
    exon_spec: Optional[float]
    exon_f1: Optional[float]
    exon_edit: Optional[float]

    utr_num_true_pos: int
    utr_num_false_pos: int
    utr_num_false_neg: int
    utr_sens: Optional[float]
    utr_spec: Optional[float]
    utr_f1: Optional[float]
    utr_edit: Optional[float]

    cds_nucl_match_coef: Optional[float]
    cds_nucl_corr_coef: Optional[float]
    cds_nucl_sens: Optional[float]
    cds_nucl_spec: Optional[float]
    cds_nucl_f1: Optional[float]
    cds_nucl_edit: Optional[float]

    utr_nucl_match_coef: Optional[float]
    utr_nucl_corr_coef: Optional[float]
    utr_nucl_sens: Optional[float]
    utr_nucl_spec: Optional[float]
    utr_nucl_f1: Optional[float]
    utr_nucl_edit: Optional[float]

    overall_match_coef: Optional[float]

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
    def from_parseval_summary(
        cls,
        summary: ParsEvalSummary
    ) -> "AnnotationStatistics":
        return cls(
            summary.total_comparisons.perfect_matches.num_matches,
            summary.total_comparisons.cds_matches.num_matches,
            summary.total_comparisons.exon_matches.num_matches,
            summary.total_comparisons.utr_matches.num_matches,
            summary.cds_structure.num_pred_match_ref,
            summary.cds_structure.num_pred_dont_match_ref,
            summary.cds_structure.num_ref_dont_match_pred,
            summary.cds_structure.sensitivity,
            summary.cds_structure.specificity,
            summary.cds_structure.f1_score,
            summary.cds_structure.edit_distance,
            summary.exon_structure.num_pred_match_ref,
            summary.exon_structure.num_pred_dont_match_ref,
            summary.exon_structure.num_ref_dont_match_pred,
            summary.exon_structure.sensitivity,
            summary.exon_structure.specificity,
            summary.exon_structure.f1_score,
            summary.exon_structure.edit_distance,
            summary.utr_structure.num_pred_match_ref,
            summary.utr_structure.num_pred_dont_match_ref,
            summary.utr_structure.num_ref_dont_match_pred,
            summary.utr_structure.sensitivity,
            summary.utr_structure.specificity,
            summary.utr_structure.f1_score,
            summary.utr_structure.edit_distance,
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
            summary.overall_nucleotide_comparison.matching_coefficient,
        )

    @classmethod
    def from_obj(cls, d: Mapping[Any, Any]) -> "AnnotationStatistics":
        out = dict()
        for column in cls.columns():
            out[column] = d[column]
        return cls(**out)

    def to_obj(self) -> Dict[str, Union[int, float]]:
        out = dict()
        for column in self.columns():
            out[column] = getattr(self, column)
        return out


class Result(NamedTuple):

    config: HintConfig
    result: Optional[AnnotationStatistics]
    errors: List[str]

    @classmethod
    def from_obj(cls, d: Mapping[Any, Any]) -> "Result":
        config = HintConfig.from_obj(d["config"])
        result = fmap(
            AnnotationStatistics.from_obj,
            d.get("result", None)
        )

        errors = d.get("errors", [])
        if not isinstance(errors, list):
            raise ValueError("Errors is expected to be a list.")

        if not all(isinstance(e, str) for e in errors):
            raise ValueError("Errors is expected to be a list of strings.")

        return cls(config, result, errors)

    def to_obj(self):
        config = self.config.to_obj()
        result = self.result.to_obj()
        return {"config": config, "result": result, "errors": self.errors}
