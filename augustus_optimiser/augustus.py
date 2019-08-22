from typing import Tuple
from typing import Union

from Bio.Application import AbstractCommandline, _Option, _Argument
from typing import Optional


def check_option_value(valid: Tuple[str], string: str) -> None:
    assert isinstance(string, str)
    if string not in valid:
        raise ValueError(
            f"Invalid option value '{string}' received. "
            f"Expected one of {repr(valid)}."
        )
    return


def check_augustus_strand(string: str) -> None:
    valid = ("both", "forward", "backward")
    check_option_value(valid, string)
    return


def check_augustus_gene_model(string: str) -> None:
    valid = ("partial", "intronless", "complete", "atleastone", "exactlyone")
    check_option_value(valid, string)
    return


def check_augustus_on_off(string: str) -> None:
    valid = ("on", "off")
    check_option_value(valid, string)
    return


def check_augustus_true_false(string: str) -> None:
    valid = ("true", "false")
    check_option_value(valid, string)
    return


def check_is_str(string: str) -> None:
    if not isinstance(string, str):
        raise ValueError(
            f"The parameter must be a string. You provided '{string}'."
        )


def check_is_int(num: int) -> None:
    if not isinstance(num, int):
        raise ValueError(
            f"The parameter must be an integer. You provided '{num}'."
        )


def check_is_int_in_range(num: int, min: int, max: int) -> None:
    check_is_int(num)
    if not (min <= num <= max):
        raise ValueError(
            f"The parameter must be between {min} and {max}. "
            f"You provided '{num}'."
        )


def check_is_numeric(num: Union[int, float]) -> None:
    if not isinstance(num, (int, float)):
        raise ValueError(
            f"The parameter must be a number. You provided '{num}'."
        )


def check_is_probability(num: Union[int, float]) -> None:
    check_is_numeric(num)
    if not (0 <= num <= 1):
        raise ValueError(
            "The parameter must be a number between 0 and 1. "
            f"You provided '{num}'."
        )


class Augustus(AbstractCommandline):

    def __init__(self, cmd="augustus", **kwargs):

        self.parameters = [
            _Option(["--species", "species"],
                    "The species to use from the config path.",
                    is_required=True),
            _Option(["--strand", "strand"],
                    ("Report predicted genes on 'both' strands, "
                     "just the 'forward', or just the 'backward' strand. "
                     "default is 'both'"),
                    checker_function=check_augustus_strand),
            _Option(["--genemodel", "genemodel"],
                    ("'partial': allow prediction of incomplete genes at the "
                     "sequence boundaries (default). 'intronless': only "
                     "predict single-exon genes like in prokaryotes and some "
                     "eukaryotes. 'complete': only predict complete genes "
                     "'atleastone': predict at least one complete gene. "
                     "'exactlyone': predict exactly one complete gene"),
                    checker_function=check_augustus_gene_model),
            _Option(["--singlestrand", "singlestrand"],
                    ("predict genes independently on each strand, allow "
                     "overlapping genes on opposite strands This option is "
                     "turned off by default."),
                    checker_function=check_augustus_true_false),
            _Option(["--hintsfile", "hintsfile"],
                    ("When this option is used the prediction considering "
                     "hints (extrinsic information) is turned on. "
                     "hintsfilename contains the hints in gff format."),
                    checker_function=check_is_str,
                    filename=True),
            _Option(["--extrinsicCfgFile", "extrinsic_cfg_file"],
                    ("This file contains the list of used sources for the "
                     "hints and their boni and mali. If not specified the "
                     "file 'extrinsic.cfg' in the config directory "
                     "$AUGUSTUS_CONFIG_PATH is used."),
                    checker_function=check_is_str,
                    filename=True),
            _Option(["--maxDNAPieceSize", "max_dna_piece_size"],
                    ("This value specifies the maximal length of the pieces "
                     "that the sequence is cut into for the core algorithm "
                     "(Viterbi) to be run. Default is "
                     "--maxDNAPieceSize=200000. AUGUSTUS tries to place the "
                     "boundaries of these pieces in the intergenic region, "
                     "which is inferred by a preliminary prediction. "
                     "GC-content dependent parameters are chosen for each "
                     "piece of DNA if /Constant/decomp_num_steps > 1 for "
                     "that species. This is why this value should not be set "
                     "very large, even if you have plenty of memory."),
                    checker_function=check_is_int),
            _Option(["--protein", "protein"],
                    "Include the proteins in the output.",
                    checker_function=check_augustus_on_off),
            _Option(["--introns", "introns"],
                    "Include the introns in the output.",
                    checker_function=check_augustus_on_off),
            _Option(["--start", "start"],
                    "Include the start in the output.",
                    checker_function=check_augustus_on_off),
            _Option(["--stop", "stop"],
                    "Include the stop in the output.",
                    checker_function=check_augustus_on_off),
            _Option(["--cds", "cds"],
                    "Include the cds in the output.",
                    checker_function=check_augustus_on_off),
            _Option(["--codingseq", "codingseq"],
                    ("Include the coding seq in the output. "
                     "Unsure how this is different to --cds."),
                    checker_function=check_augustus_on_off),
            _Option(["--stopCodonExcludedFromCDS",
                     "stop_codon_excluded_from_cds"],
                    "Exclude the stop codon from cds output.",
                    checker_function=check_augustus_true_false),
            _Option(["--AUGUSTUS_CONFIG_PATH", "augustus_config_path"],
                    ("Path to config directory (if not specified as "
                     "environment variable)"),
                    checker_function=check_is_str,
                    filename=True),
            _Option(["--alternatives-from-evidence",
                     "alternatives_from_evidence"],
                    ("Report alternative transcripts when they are "
                     "suggested by hints."),
                    checker_function=check_augustus_true_false),
            _Option(["--alternatives-from-sampling",
                     "alternatives_from_sampling"],
                    ("Report alternative transcripts generated through "
                     "probabilistic sampling."),
                    checker_function=check_augustus_true_false),
            _Option(["--sample", "sample"],
                    ("The posterior probabilities are estimated using a "
                     "sampling algorithm. The parameter `--sample=n` adjusts "
                     "the number of sampling iterations. The higher 'n' is "
                     "the more accurate is the estimation but it usually "
                     "isn't important that the posterior probability is "
                     "very accurate. Every 30 sample iterations take about "
                     "the same time as one run without sampling, e.g. "
                     "`--sample=60` takes about 3 times as much time as "
                     "`--sample=0`. The default is `--sample=100`"),
                    checker_function=check_is_int),
            _Option(["--minexonintronprob", "minexonintronprob"],
                    ("When deriving alternatives_from_sampling. The posterior "
                     "probability of every exon and every intron in a "
                     "transcript must be at least 'minexonintronprob', "
                     "otherwise the transcript is not reported. "
                     "minexonintronprob=0.1 is a reasonable value."),
                    checker_function=check_is_probability),
            _Option(["--minmeanexonintronprob", "minmeanexonintronprob"],
                    ("When deriving alternatives_from_sampling. The "
                     "geometric mean of the probabilities of all exons "
                     "and introns must be at least 'minmeanexonintronprob'. "
                     "minmeanexonintronprob=0.4 is a reasonable value."),
                    checker_function=check_is_probability),
            _Option(["--maxtracks", "maxtracks"],
                    ("The maximum number of tracks when displayed in a "
                     "genome browser is 'maxtracks' (unless maxtracks=-1, "
                     "then it is unbounded). In cases where all transcripts "
                     "of a gene overlap at some position this is also the "
                     "maximal number of transcripts for that gene. I "
                     "recommend increasing the parameter 'maxtracks' for "
                     "improving sensitivity and setting 'maxtracks' to 1 "
                     "and increasing minmeanexonintronprob and/or "
                     "minexonintronprob in order to improve specificity."),
                    checker_function=check_is_int),
            _Option(["--temperature", "temperature"],
                    ("The probabilistic model of AUGUSTUS can be seen as a "
                     "rough approximation to reality. A consequence is that "
                     "the posterior probabilities for the strong exons "
                     "(e.g. the ones called by the Viterbi algorithm) tend "
                     "to be larger than the actually measured precision "
                     "(specificity) values. For example, in human, "
                     "only 94.5% of the exons with predicted posterior "
                     "probability >= 98% (under the default --sample=100) "
                     "are actually correct. If the aim of the "
                     "sampling is to produce a diverse, sensitive set of "
                     "gene structures, you can use the `--temperature=t` "
                     "parameter. Where t is one of 0,1,2,3,4,5,6,7. All "
                     "probabilities of the model are then taken to the "
                     "power of (8-t)/8, i.e. t=0 (the default) does nothing. "
                     "The larger t the more alternatives are sampled. "
                     "t=3 is a good compromise between getting a high "
                     "sensitivity but not getting too many exons sampled "
                     "in total. For t=3, 96.1% of human exons with "
                     "posterior probability >= 98% are correct."),
                    checker_function=lambda x: check_is_int_in_range(x, 0, 7)),
        ]
        proteinprofile: Optional[str] = None
        outfile: Optional[str] = None
        prediction_start: Optional[int] = None
        prediction_end: Optional[int] = None
        no_inframe_stop: bool = False
        noprediction: bool = False
        unique_gene_id: bool = False
        max_dna_piece_size: int = 200000
        progress: bool = False
        gff3: AugustusOnOff = AugustusOnOff.ON
        utr: AugustusOnOff = AugustusOnOff.OFF
        introns: AugustusOnOff = AugustusOnOff.ON
        start: AugustusOnOff = AugustusOnOff.ON
        stop: AugustusOnOff = AugustusOnOff.ON
        protein: AugustusOnOff = AugustusOnOff.OFF
        cds: AugustusOnOff = AugustusOnOff.OFF
        stop_codon_excluded_from_cds: bool = True
        codingseq: AugustusOnOff = AugustusOnOff.OFF
        contentmodels: AugustusOnOff = AugustusOnOff.ON
        return
