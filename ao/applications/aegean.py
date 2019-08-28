from Bio.Application import AbstractCommandline
from Bio.Application import _Option, _Argument, _Switch


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


class ParsEval(AbstractCommandline):

    def __init__(self, cmd="parseval", **kwargs):
        """ Construct and evaluate aegean parseval commands.

        Example
        >>> x = ParsEval(outformat="text",
        ...              reference="ref.gff3",
        ...              prediction="pred.gff3")
        >>> print(x)
        parseval --outformat text ref.gff3 pred.gff3
        """

        self.program_name = cmd
        self.parameters = [
            _Switch(["--debug", "debug"],
                    "Print debugging messages"),
            _Switch(["--help", "help"],
                    "Print help message and exit"),
            _Option(["--delta", "delta"],
                    "Extend gene loci by this many nucleotides.",
                    checker_function=check_is_positive_int,
                    equate=False,
                    ),
            _Switch(["--verbose", "verbose"],
                    "Print verbose warning messages"),
            _Switch(["--version", "version"],
                    "Print version number and exit"),
            _Option(["--datashare", "datashare"],
                    ("Indicate desired output format; possible "
                     "options: 'csv', 'text', or 'html' "
                     "(default='text'); in 'text' or 'csv' mode, "
                     "will create a single file; in 'html' mode, "
                     "will create a directory."),
                    checker_function=check_is_str,
                    filename=True,
                    equate=False,
                    ),
            _Option(["--outformat", "outformat"],
                    ("Indicate desired output format; possible "
                     "options: 'csv', 'text', or 'html' "
                     "(default='text'); in 'text' or 'csv' mode, "
                     "will create a single file; in 'html' mode, "
                     "will create a directory."),
                    checker_function=check_is_str,
                    filename=True,
                    equate=False,
                    ),
            _Switch(["--nogff3", "nogff3"],
                    ("Do no print GFF3 output corresponding to each "
                     "comparison.")),
            _Option(["--outfile", "outfile"],
                    ("File/directory to which output will be "
                     "written; default is the terminal (STDOUT)."),
                    checker_function=check_is_str,
                    filename=True,
                    equate=False,
                    ),
            _Switch(["--nopng", "nopng"],
                    ("In HTML output mode, skip generation of PNG "
                     "graphics for each gene locus.")),
            _Switch(["--summary", "summary"],
                    ("Only print summary statistics, do not print "
                     "individual comparisons.")),
            _Switch(["--overwrite", "overwrite"],
                    "Force overwrite of any existing output files."),
            _Option(["--refrlabel", "refrlabel"],
                    "Optional label for reference annotations.",
                    checker_function=check_is_str,
                    equate=False,
                    ),
            _Option(["--predlabel", "predlabel"],
                    "Optional label for prediction annotations.",
                    checker_function=check_is_str,
                    equate=False,
                    ),
            _Switch(["--makefilter", "makefilter"],
                    ("Create a default configuration file for "
                     "filtering reported results and quit, "
                     "performing no comparisons.")),
            _Option(["--filterfile", "filterfile"],
                    ("Use the indicated configuration file to "
                     "filter reported results."),
                    checker_function=check_is_str,
                    equate=False,
                    ),
            _Option(["--maxtrans", "maxtrans"],
                    ("Maximum transcripts allowed per locus; use 0 "
                     "to disable limit; default is 32"),
                    checker_function=check_is_positive_int,
                    equate=False,
                    ),
            _Argument(["reference"],
                      "The reference gff3.",
                      checker_function=check_is_str,
                      filename=True,
                      is_required=True),
            _Argument(["prediction"],
                      "The prediction gff3.",
                      checker_function=check_is_str,
                      filename=True,
                      is_required=True),
        ]

        super().__init__(cmd, **kwargs)
        return
