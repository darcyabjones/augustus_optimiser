import yaml
import json
import argparse

from tempfile import NamedTemporaryFile

from typing import Optional
from typing import Set, List, Sequence

from functools import partial

from mpi4py import MPI
from mpi4py.futures import MPIPoolExecutor

from ao.parsers import extrinsic_metaparameters
from ao.stats import AnnotationStatistics, Result
from ao.hints import HintConfig, HintConfigFactory
from ao.applications.augustus import Augustus
from ao.applications.genometools import GFF3, Merge
from ao.applications.aegean import ParsEval
from ao.parsers.parseval_summary import ParsEvalSummary


EXTRINSIC_DESC = "Generate and evaluate many extrinsic hint parameters."


def cli_extrinsic(parser):
    parser.add_argument(
        "config",
        type=argparse.FileType('r'),
        help="Input config file."
    )

    parser.add_argument(
        "genome",
        type=str,
        help="Input genome fasta file."
    )

    parser.add_argument(
        "gff",
        type=str,
        help="Input 'true' gff file."
    )

    parser.add_argument(
        "hints",
        type=str,
        help="Input hints file."
    )

    parser.add_argument(
        "-f", "--config-format",
        dest="config_format",
        choices=["json", "yaml"],
        type=str,
        default="yaml",
        help="What format is the config file in?"
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType("w"),
        help="Output line delimited json results file."
    )

    parser.add_argument(
        "-t", "--tmpdir",
        type=str,
        default=None,
        help="Where to write temporary files."
    )

    parser.add_argument(
        "-c", "--continue",
        dest="continue_",
        type=argparse.FileType("r"),
        help="Continue using an existing result file."
    )

    parser.add_argument(
        "-n", "--ntrials",
        type=int,
        default=1000,
        help="How many random configs should be tried"
    )

    parser.add_argument(
        "--species",
        type=str,
        required=True,
        help=(
            "The species id to use for augustus. "
            "Passed to the augustus `--species` parameter."
        )
    )

    parser.add_argument(
        "--no-softmasking",
        dest="softmasking",
        action="store_false",
        default=True,
        help=("Prevent augustus from interpreting lowercase bases as "
              "softmasked. Sets augustus parameter `--softmasking=off`.")
    )

    parser.add_argument(
        "--utr",
        action="store_true",
        default=False,
        help=("Predict UTRs using augustus. Requires that the UTR "
              "model is trained for the specified species. "
              "Sets augustus parameter `--UTR=on`.")
    )

    parser.add_argument(
        "--alt-from-evidence",
        action="store_true",
        default=False,
        help=("Tell augustus to predict alternative transcripts using "
              "evidence. Sets augustus parameter "
              "`--alternatives-from-evidence=true`.")
    )

    parser.add_argument(
        "--min-intron-len",
        type=int,
        default=None,
        help=("Set the minimum intron length that augustus should predict. "
              "Sets the augustus parameter `--min_intron_len`.")
    )

    parser.add_argument(
        "--allow-hinted-splicesites",
        type=str,
        default=None,
        help=("Allow non-standard splicesites when supported by hints. "
              "Provide as a comma separated list. e.g. gcag,atac. "
              "Passed to augustus parameter `--allow_hinted_splicesites`.")
    )

    parser.add_argument(
        "--singlestrand",
        action="store_true",
        default=False,
        help=("Tell augustus to look for overlapping genes in "
              "reverse orientation. For non-utr runs, this passes the "
              "augustus parameter `--singlestrand`. If the `--utr` option is "
              "used augustus would fail because it is not implemented, "
              "instead we predict each strand individually and merge the "
              "results.")
    )

    parser.add_argument(
        "--augustus-config-path",
        type=str,
        default=None,
        help=("Path to the augustus config directory. By default will be "
              "determined from the AUGUSTUS_CONFIG_PATH environment variable "
              "by augustus. Passed to `--AUGUSTUS_CONFIG_PATH` parameter ")
    )

    # stop_codon_excluded_from_cds="false",

    return


def parse_existing_trials(handle: Sequence[str]) -> List[Result]:
    results = list()
    for line in handle:
        d = json.loads(line)
        result = Result.from_obj(d)
        results.append(result)

    return results


def generate_configs(
    n: int,
    seen: Set[str],
    factory: HintConfigFactory
) -> List[HintConfig]:

    out: List[HintConfig] = list()

    while len(out) < n:
        cfg = factory.get()
        str_cfg = str(cfg)
        if str_cfg not in seen:
            seen.add(str_cfg)
            out.append(cfg)

    return out


def write_config(config: HintConfig, tmpdir: Optional[str]) -> str:
    config_file_name = write_to_named_tmp(str(config), "config_", tmpdir)
    return config_file_name


def write_to_named_tmp(string: str, prefix: str, tmpdir: Optional[str]) -> str:
    try:
        tmpfile = NamedTemporaryFile(
            mode="w",
            prefix=prefix,
            delete=False,
            dir=tmpdir
        )

        print(string, file=tmpfile)

        tmpfile_name = tmpfile.name
    finally:
        tmpfile.close()

    return tmpfile_name


def filter_augustus_gff(lines: Sequence[str]) -> str:
    out = list()
    for line in lines:
        if line.startswith("#"):
            continue

        split_line = line.split("\t")

        if len(split_line) < 3:
            continue

        if split_line[2] == "transcript":
            split_line[2] = "mRNA"
            out_line = "\t".join(split_line)
        else:
            out_line = line

        out.append(out_line)

    return '\n'.join(out)


def run_augustus(
    genome: str,
    species: str,
    config: str,
    hints: str,
    softmasking: bool,
    utr: bool,
    alt_from_evidence: bool,
    min_intron_len: Optional[int],
    allow_hinted_splicesites: Optional[str],
    singlestrand: bool,
    augustus_config_path: Optional[str],
    tmpdir: Optional[str],
    strand: str = "both",
) -> List[str]:
    """ """

    if strand not in ("forward", "backward", "both"):
        raise ValueError("Invalid augustus strand.")

    if singlestrand and utr:
        run_augustus_strand = partial(
            run_augustus,
            genome=genome,
            species=species,
            config=config,
            hints=hints,
            softmasking=softmasking,
            utr=utr,
            alt_from_evidence=alt_from_evidence,
            min_intron_len=min_intron_len,
            allow_hinted_splicesites=allow_hinted_splicesites,
            singlestrand=False,
            tmpdir=tmpdir,
            augustus_config_path=augustus_config_path,
        )
        fwd = run_augustus_strand(strand="forward")
        rev = run_augustus_strand(strand="backward")
        return fwd + rev

    softmask_onoff = "on" if softmasking else "off"
    utr_onoff = "on" if utr else "off"
    alt_from_evidence_truefalse = "true" if alt_from_evidence else "false"
    singlestrand_truefalse = "true" if singlestrand else "false"

    augustus_cmd = Augustus(
        infile=genome,
        strand=strand,
        species=species,
        extrinsic_cfg_file=config,
        hintsfile=hints,
        softmasking=softmask_onoff,
        singlestrand=singlestrand_truefalse,
        utr=utr_onoff,
        alternatives_from_evidence=alt_from_evidence_truefalse,
        gff3="on",
        progress="false",
        protein="off",
        start="off",
        stop="off",
        introns="off",
        stop_codon_excluded_from_cds="false",
    )

    if min_intron_len is not None:
        augustus_cmd.min_intron_len = min_intron_len

    if allow_hinted_splicesites is not None:
        augustus_cmd.allow_hinted_splicesites = allow_hinted_splicesites

    if augustus_config_path is not None:
        augustus_cmd.augustus_config_path = augustus_config_path

    stdout, stderr = augustus_cmd(stdout=True, stderr=True)
    # Also need to be able to catch errors, though we don't seem to be
    # encountering them really.

    filtered = filter_augustus_gff(stdout.split('\n'))
    file_name = write_to_named_tmp(filtered, "augustus_", tmpdir)
    return [file_name]


def evaluate(
    config: HintConfig,
    genome: str,
    gff: str,
    hints: str,
    species: str,
    tmpdir: Optional[str],
    singlestrand: bool,
    softmasking: bool,
    utr: bool,
    alt_from_evidence: bool,
    min_intron_len: Optional[int],
    allow_hinted_splicesites: Optional[str],
    augustus_config_path: Optional[str],
) -> Result:

    config_file_name = write_config(config, tmpdir)

    augustus_gffs = run_augustus(
        genome,
        species,
        config_file_name,
        hints,
        softmasking,
        utr,
        alt_from_evidence,
        min_intron_len,
        allow_hinted_splicesites,
        singlestrand,
        augustus_config_path,
        tmpdir,
        strand="both",
    )

    tidied_gffs = []
    for augustus_gff in augustus_gffs:
        tidied_gff = augustus_gff + "_genometools"
        genometools_cmd = GFF3(
            infile=augustus_gff,
            outfile=tidied_gff,
            sort=True,
            tidy=True,
            addintrons=False,
        )

        gt_stdout, gt_stderr = genometools_cmd()
        tidied_gffs.append(tidied_gff)

    if len(tidied_gffs) == 1:
        genometools_outfile = tidied_gffs[0]
    else:
        genometools_outfile = tidied_gffs[0] + "_merged"

        gtmerge_cmd = Merge(
            infiles=tidied_gffs,
            outfile=genometools_outfile,
            tidy=True,
        )
        gtm_stdout, gtm_stderr = gtmerge_cmd()

    aegean_cmd = ParsEval(
        outformat="text",
        nogff3=True,
        reference=gff,
        summary=True,
        prediction=genometools_outfile,
    )

    aegean_stdout, aegean_stderr = aegean_cmd(stdout=True, stderr=True)

    # Iter is necessary to make the list behave like an iterator.
    # ie. it gets consumed.
    parsed_result = ParsEvalSummary.parse(iter(aegean_stdout.split('\n')))
    ann_stats = AnnotationStatistics.from_parseval_summary(parsed_result)

    return Result(config, ann_stats, [])


def run_extrinsic(args) -> None:

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # The pool executor handles all other workers.
    if rank != 0:
        return

    if args.config_format == "yaml":
        cfg = yaml.safe_load(args.config)
    elif args.config_format == "json":
        cfg = json.load(args.config)
    else:
        raise ValueError("Currently only json and yaml formats are supported.")

    cfg_factory = extrinsic_metaparameters.parser(cfg)

    if args.continue_ is None:
        results: List[Result] = []
    else:
        results = parse_existing_trials(args.continue_)

    seen = {str(result.config) for result in results}

    configs = generate_configs(args.ntrials, seen, cfg_factory)

    run = partial(
        evaluate,
        genome=args.genome,
        gff=args.gff,
        hints=args.hints,
        species=args.species,
        tmpdir=args.tmpdir,
        softmasking=args.softmasking,
        singlestrand=args.singlestrand,
        utr=args.utr,
        alt_from_evidence=args.alt_from_evidence,
        min_intron_len=args.min_intron_len,
        allow_hinted_splicesites=args.allow_hinted_splicesites,
        augustus_config_path=args.augustus_config_path,
    )

    # max_workers=None means that the universe size is defined by the master
    # mpi process.
    executor = MPIPoolExecutor(max_workers=None)

    for result in executor.map(run, configs, unordered=True):
        results.append(result)
        seen.add(str(result.config))
        print(json.dumps(result.to_obj()), file=args.outfile)

    return
