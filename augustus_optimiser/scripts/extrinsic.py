import yaml
import json
import argparse

from tempfile import NamedTemporaryFile

from typing import Optional
from typing import Set, List, Sequence

from augustus_optimiser.parsers import extrinsic_metaparameters
from augustus_optimiser.stats import AnnotationStatistics, Result
from augustus_optimiser.hints import HintConfig, HintConfigFactory
from augustus_optimiser.applications.augustus import Augustus
from augustus_optimiser.applications.genometools import GFF3
from augustus_optimiser.applications.aegean import ParsEval
from augustus_optimiser.parsers.parseval_summary import ParsEvalSummary


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
        help="The species to use for augustus."
    )

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
    try:
        config_file = NamedTemporaryFile(
            mode="w",
            prefix="config_",
            delete=False,
            dir=tmpdir
        )

        print(str(config), file=config_file)

        config_file_name = config_file.name
    finally:
        config_file.close()

    return config_file_name


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


def run_augustus(genome, species, config, hints) -> str:
    """ """

    augustus_cmd = Augustus(
        infile=genome,
        species=species,
        extrinsic_cfg_file=config,
        hintsfile=hints,
        softmasking="on",
        gff3="on",
        utr="off",
        progress="false",
        allow_hinted_splicesites="gcag,atac,ctac,gaag",
        alternatives_from_evidence="true",
        protein="off",
        min_intron_len=20,
        stop_codon_excluded_from_cds="false",
    )

    stdout, stderr = augustus_cmd(stdout=True, stderr=True)
    # Also need to be able to catch errors.

    filtered = filter_augustus_gff(stdout.split('\n'))
    return filtered


def evaluate(
    genome: str,
    gff: str,
    hints: str,
    species: str,
    config: HintConfig,
    tmpdir: Optional[str]
) -> Result:

    config_file_name = write_config(config, tmpdir)

    augustus_gff = run_augustus(genome, species, config_file_name, hints)

    genometools_outfile = config_file_name + "_genometools"
    genometools_cmd = GFF3(
        infile="-",
        outfile=genometools_outfile,
        sort=True,
        tidy=True,
        addintrons=True,
    )

    gt_stdout, gt_stderr = genometools_cmd(stdin=augustus_gff)

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


def extrinsic(args) -> None:

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

    for config in configs:
        result = evaluate(
            args.genome,
            args.gff,
            args.hints,
            args.species,
            config,
            args.tmpdir
        )
        print(json.dumps(result.to_obj()), file=args.outfile)

    # scatter_options
    # Write out finished as you get them.
    return
