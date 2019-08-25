from typing import List
from typing import Union

from Bio.Application import AbstractCommandline, _Option, _Argument


def check_option_value(valid: List[str], string: str) -> bool:
    assert isinstance(string, str)
    if string not in valid:
        raise ValueError(
            f"Invalid option value '{string}' received. "
            f"Expected one of {repr(valid)}."
        )
    return True


def check_augustus_strand(string: str) -> bool:
    valid = ["both", "forward", "backward"]
    check_option_value(valid, string)
    return True


def check_augustus_gene_model(string: str) -> bool:
    valid = ["partial", "intronless", "complete", "atleastone", "exactlyone"]
    check_option_value(valid, string)
    return True


def check_augustus_on_off(string: str) -> bool:
    valid = ["on", "off"]
    check_option_value(valid, string)
    return True


def check_augustus_true_false(string: str) -> bool:
    valid = ["true", "false"]
    check_option_value(valid, string)
    return True


def check_is_str(string: str) -> bool:
    if not isinstance(string, str):
        raise ValueError(
            f"The parameter must be a string. You provided '{string}'."
        )
    return True


def check_is_int(num: int) -> bool:
    if not isinstance(num, int):
        raise ValueError(
            f"The parameter must be an integer. You provided '{num}'."
        )
    return True


def check_is_positive_int(num: int) -> bool:
    check_is_int(num)
    if num < 0:
        raise ValueError(
            f"The parameter cannot be negative. You provided '{num}'."
        )
    return True


def check_is_int_in_range(num: int, min: int, max: int) -> bool:
    check_is_int(num)
    if not (min <= num <= max):
        raise ValueError(
            f"The parameter must be between '{min}' and '{max}'. "
            f"You provided '{num}'."
        )
    return True


def check_is_numeric(num: Union[int, float]) -> bool:
    if not isinstance(num, (int, float)):
        raise ValueError(
            f"The parameter must be a number. You provided '{num}'."
        )
    return True


def check_is_positive_numeric(num: Union[int, float]) -> bool:
    check_is_numeric(num)
    if num < 0:
        raise ValueError(
            f"The parameter cannot be negative. You provided '{num}'."
        )
    return True


def check_is_numeric_in_range(
    num: Union[int, float],
    min: Union[int, float],
    max: Union[int, float]
) -> bool:
    check_is_numeric(num)
    if not (min <= num <= max):
        raise ValueError(
            f"The parameter must be between '{min}' and '{max}'. "
            f"You provided '{num}'."
        )
    return True


def check_is_probability(num: Union[int, float]) -> bool:
    check_is_numeric_in_range(num, 0, 1)
    return True


class Augustus(AbstractCommandline):

    def __init__(self, cmd="augustus", **kwargs):
        self.program_name = cmd
        self.parameters = [
            _Argument(["infile"],
                      "The input fasta or genbank file to run analysis on.",
                      checker_function=check_is_str,
                      filename=True,
                      is_required=True),
            _Option(["--species", "species"],
                    "The species to use from the config path.",
                    is_required=True,
                    checker_function=check_is_str),
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
                    checker_function=check_is_positive_int),
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
                    checker_function=check_is_positive_int),
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
            _Option(["--proteinprofile", "proteinprofile"],
                    "Read a protein profile from file filename.",
                    checker_function=check_is_str,
                    filename=True),
            _Option(["--predictionStart", "prediction_start"],
                    ("Define the range of the sequence for which predictions "
                     "should be found."),
                    checker_function=check_is_positive_int),
            _Option(["--predictionEnd", "prediction_end"],
                    ("Define the range of the sequence for which predictions "
                     "should be found."),
                    checker_function=check_is_positive_int),
            _Option(["--gff3", "gff3"],
                    "Output in gff3 format",
                    checker_function=check_augustus_on_off),
            _Option(["--UTR", "utr"],
                    ("Predict the untranslated regions in addition to "
                     "the coding sequence."),
                    checker_function=check_augustus_on_off),
            _Option(["--outfile", "outfile"],
                    "print output to filename instead to standard output.",
                    checker_function=check_is_str,
                    filename=True),
            _Option(["--noInFrameStop", "no_in_frame_stop"],
                    ("Don't report transcripts with in-frame stop codons. "
                     "Otherwise, intron-spanning stop codons could occur."),
                    checker_function=check_augustus_true_false),
            _Option(["--noprediction", "noprediction"],
                    ("If true and input is in genbank format, no "
                     "prediction is made. Useful for getting the annotated "
                     "protein sequences."),
                    checker_function=check_augustus_true_false),
            _Option(["--contentmodels", "contentmodels"],
                    ("If 'off' the content models are disabled (all "
                     "emissions uniformly 1/4). The content models are; "
                     "coding region Markov chain (emiprobs), initial k-mers "
                     "in coding region (Pls), intron and intergenic regin "
                     "Markov chain. This option is intended for special "
                     "applications that require judging gene structures "
                     "from the signal models only, e.g. for predicting "
                     "the effect of SNPs or mutations on splicing. For "
                     "all typical gene predictions, this should be true. "
                     "Default: on"),
                    checker_function=check_augustus_on_off),
            _Option(["--progress", "progress"],
                    "Print running progress reports.",
                    checker_function=check_augustus_true_false),
            _Option(["--uniqueGeneId", "unique_gene_id"],
                    "If true, output gene identifyers like this: seqname.gN",
                    checker_function=check_augustus_true_false),
            _Option(["--allow_hinted_splicesites", "allow_hinted_splicesites"],
                    ("Allow non-standard splice sites when supported by "
                     "hints. Provide as a comma separated list e.g. atac"),
                    checker_function=check_is_str),
            _Option(["--alnfile", "alnfile"],
                    "The multiple genome alignment file to use for CPG.",
                    checker_function=check_is_str,
                    filename=True),
            _Option(["--min_intron_len", "min_intron_len"],
                    "The minimum allowable intron length in bp.",
                    checker_function=check_is_positive_int),
            _Option(["--softmasking", "softmasking"],
                    "Softmask sequences.",
                    checker_function=check_augustus_on_off),
            _Option(["--errfile", "errfile"],
                    "Output stderr to a file.",
                    checker_function=check_is_str,
                    filename=True),

        ]
        super().__init__(cmd, **kwargs)
        return


"""
Remaining augustus options.
Unsure what they do, some are only used for training, some for CPG mode.

"/augustus/verbosity",
"/BaseCount/weighingType",
"/BaseCount/weightMatrixFile",
"bridge_genicpart_bonus",
"canCauseAltSplice",
"capthresh",
"checkExAcc",
"codonAlignmentFile",
"complete_genes",
"/CompPred/assmotifqthresh",
"/CompPred/assqthresh",
"/CompPred/dssqthresh",
"/CompPred/compSigScoring",
"/CompPred/conservation",
"/CompPred/covPen",
"/CompPred/ec_score",
"/CompPred/exon_gain",
"/CompPred/exon_loss",
"/CompPred/ali_error",
"/CompPred/computeNumSubs",
"/CompPred/maxIterations",
"/CompPred/mil_factor",
"/CompPred/ec_factor",
"/CompPred/ec_addend",
"/CompPred/ec_thold",
"/CompPred/ic_thold",
"/CompPred/genesWithoutUTRs",
"/CompPred/liftover_all_ECs",
"/CompPred/logreg",
"/CompPred/maxCov",
"/CompPred/omega",
"/CompPred/num_omega",
"/CompPred/num_features",
"/CompPred/scale_codontree",
"/CompPred/oeExtensionWidth",
"/CompPred/onlySampledECs",
"/CompPred/only_species",
"/CompPred/outdir_orthoexons",
"/CompPred/outdir",
"/CompPred/printOrthoExonAli",
"/CompPred/printConservationWig",
"/CompPred/phylo_factor",
"/CompPred/phylo_model",
"/CompPred/dd_factor",
"/CompPred/dd_rounds",
"/CompPred/dd_step_rule",
"/CompPred/dualdecomp",
"/CompPred/overlapcomp",
"/CompPred/lambda",
"/Constant/almost_identical_maxdiff",
"/Constant/amberprob",
"/Constant/ass_end",
"/Constant/ass_maxbinsize",
"/Constant/ass_start",
"/Constant/ass_upwindow_size",
"/Constant/decomp_num_at",
"/Constant/decomp_num_gc",
"/Constant/decomp_num_steps",
"/Constant/dss_end",
"/Constant/dss_maxbinsize",
"/Constant/dss_start",
"/Constant/gc_range_max",
"/Constant/gc_range_min",
"/Constant/init_coding_len",
"/Constant/intterm_coding_len",
"/Constant/max_contra_supp_ratio",
"/Constant/min_coding_len",
"/Constant/ochreprob",
"/Constant/opalprob",
"/Constant/probNinCoding",
"/Constant/subopt_transcript_threshold",
"/Constant/tis_maxbinsize",
"/Constant/trans_init_window",
"/Constant/tss_upwindow_size",
"CRF",
"CRF_N",
"CRF_TRAIN",
"CRFtrainCDS",
"CRFtrainIntron",
"CRFtrainIgenic",
"CRFtrainSS",
"CRFtrainUTR",
"CRFtrainTIS",
"dbaccess",
"dbhints",
"/EHMMTraining/state00",
"/EHMMTraining/state01",
"/EHMMTraining/state02",
"/EHMMTraining/state03",
"/EHMMTraining/statecount",
"/EHMMTraining/trainfile",
"emiprobs",
"evalset",
"exoncands",
"/ExonModel/etorder",
"/ExonModel/etpseudocount",
"/ExonModel/exonlengthD",
"/ExonModel/infile",
"/ExonModel/k",
"/ExonModel/lenboostE",
"/ExonModel/lenboostL",
"/ExonModel/maxexonlength",
"/ExonModel/minexonlength",
"/ExonModel/minPatSum",
"/ExonModel/minwindowcount",
"/ExonModel/outfile",
"/ExonModel/patpseudocount",
"/ExonModel/slope_of_bandwidth",
"/ExonModel/tisalpha",
"/ExonModel/tis_motif_memory",
"/ExonModel/tis_motif_radius",
"/ExonModel/verbosity",
"exonnames",
"GCwinsize",
"GD_stepsize",
"/genbank/verbosity",
"/HMMTraining/savefile",
"/IGenicModel/infile",
"/IGenicModel/k",
"/IGenicModel/outfile",
"/IGenicModel/patpseudocount",
"/IGenicModel/verbosity",
"/IntronModel/allow_dss_consensus_gc",
"/IntronModel/ass_motif_memory",
"/IntronModel/ass_motif_radius",
"/IntronModel/asspseudocount",
"/IntronModel/d",
"/IntronModel/dssneighborfactor",
"/IntronModel/dsspseudocount",
"/IntronModel/infile",
"/IntronModel/k",
"/IntronModel/minwindowcount",
"/IntronModel/non_ag_ass_prob",
"/IntronModel/non_gt_dss_prob",
"/IntronModel/outfile",
"/IntronModel/patpseudocount",
"/IntronModel/sf_with_motif",
"/IntronModel/slope_of_bandwidth",
"/IntronModel/splicefile",
"/IntronModel/ssalpha",
"/IntronModel/verbosity",
"keep_viterbi",
"label_flip_prob",
"learning_rate",
"lossweight",
"locustree",
"maxDNAPieceSize",
"maxOvlp",
"maxtracks",
"max_sgd_inter",
"mea",
"mea_evaluation",
"/MeaPrediction/no_compatible_edges",
"/MeaPrediction/alpha_E",
"/MeaPrediction/alpha_I",
"/MeaPrediction/x0_E",
"/MeaPrediction/x0_I",
"/MeaPrediction/x1_E",
"/MeaPrediction/x1_I",
"/MeaPrediction/y0_E",
"/MeaPrediction/y0_I",
"/MeaPrediction/i1_E",
"/MeaPrediction/i1_I",
"/MeaPrediction/i2_E",
"/MeaPrediction/i2_I",
"/MeaPrediction/j1_E",
"/MeaPrediction/j1_I",
"/MeaPrediction/j2_E",
"/MeaPrediction/j2_I",
"/MeaPrediction/weight_base",
"/MeaPrediction/r_be",
"/MeaPrediction/r_bi",
"/MeaPrediction/weight_exon",
"/MeaPrediction/weight_gene",
"/MeaPrediction/weight_utr",
"minexonintronprob",
"minmeanexonintronprob",
"optCfgFile",
"printHints",
"printMEA",
"printOEs",
"printSampled",
"printGeneRangesBED",
"printGeneRangesGFF",
"param_outfile",
"print_blocks",
"print_utr",
"/ProteinModel/allow_truncated",
"/ProteinModel/exhaustive_substates",
"/ProteinModel/block_threshold_spec",
"/ProteinModel/block_threshold_sens",
"/ProteinModel/blockpart_threshold_spec",
"/ProteinModel/blockpart_threshold_sens",
"/ProteinModel/global_factor_threshold",
"/ProteinModel/absolute_malus_threshold",
"/ProteinModel/invalid_score",
"/ProteinModel/weight",
"referenceFile",
"refSpecies",
"rescaleBoni",
"rLogReg",
"scorediffweight", // temp
"speciesfilenames",
"tieIgenicIntron",
"trainFeatureFile",
"translation_table",
"treefile",
"truncateMaskedUTRs",
"tss",
"tts",
"uniqueCDS",
"useAminoAcidRates",
"useNonCodingModel",
"use_sgd",
"useTFprob",
"/UtrModel/d_polya_cleavage_max",
"/UtrModel/d_polya_cleavage_min",
"/UtrModel/d_polyasig_cleavage",
"/UtrModel/d_tss_tata_max",
"/UtrModel/d_tss_tata_min",
"/UtrModel/exonlengthD",
"/UtrModel/infile",
"/UtrModel/k",
"/UtrModel/max3singlelength",
"/UtrModel/max3termlength",
"/UtrModel/maxexonlength",
"/UtrModel/minwindowcount",
"/UtrModel/outfile",
"/UtrModel/patpseudocount",
"/UtrModel/prob_polya",
"/UtrModel/slope_of_bandwidth",
"/UtrModel/tata_end",
"/UtrModel/tata_pseudocount",
"/UtrModel/tata_start",
"/UtrModel/tss_end",
"/UtrModel/tss_start",
"/UtrModel/tssup_k",
"/UtrModel/tssup_patpseudocount",
"/UtrModel/tts_motif_memory",
"/UtrModel/utr3patternweight",
"/UtrModel/utr5patternweight",
"/UtrModel/utr3prepatternweight",
"/UtrModel/utr5prepatternweight",
"/UtrModel/verbosity"
"""
