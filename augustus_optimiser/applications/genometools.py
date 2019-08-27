from Bio.Application import AbstractCommandline
from Bio.Application import _Option, _Argument, _StaticArgument, _Switch


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


class GFF3(AbstractCommandline):

    def __init__(self, cmd="gt", **kwargs):
        """ Construct and evaluate genometools gff3 commands.

        Example
        >>> x = GFF3(infile="test.gff3", sort=True)
        >>> print(x)
        gt gff3 -sort test.gff3
        """

        self.program_name = f"{cmd} gff3"
        self.parameters = [
            _StaticArgument("gff3"),
            _Switch(["-sort", "sort"],
                    ("Sort the GFF3 features (memory consumption is "
                     "proportional to the input file size(s))")
                    ),
            _Switch(["-sortlines", "sortlines"],
                    ("sort the GFF3 features on a strict line basis "
                     "(not sorted asdefined by GenomeTools).")
                    ),
            _Switch(["-sortnum", "sortnum"],
                    ("enable natural numeric sorting for sequence regions "
                     "(not sorted as defined by GenomeTools)")
                    ),
            _Switch(["-tidy", "tidy"],
                    ("Try to tidy the GFF3 files up during parsing.")
                    ),
            _Switch(["-retainids", "retainids"],
                    ("when available, use the original IDs provided in the "
                     "source file")
                    ),
            _Switch(["-checkids", "checkids"],
                    ("make sure the ID attributes are unique within the "
                     "scope of each GFF3_file, as required by GFF3 "
                     "specification. If features with the same Parent "
                     "attribute are not separated by a # line the GFF3 "
                     "parser tries to treat them as a multi-line feature. "
                     "This requires at least matching sequence IDs and types")
                    ),
            _Switch(["-addids", "dont_addids"],
                    "add missing '##sequence-region' lines automatically"
                    ),
            _Switch(["-fixregionboundaries", "fixregionboundaries"],
                    ('automatically adjust "##sequence-region" lines to '
                     "contain all their features")
                    ),
            _Switch(["-addintrons", "addintrons"],
                    "add intron features between existing exon features"
                    ),
            _Option(["-offset", "offset"],
                    "transform all features by the given offset",
                    checker_function=check_is_positive_int,
                    equate=False,
                    ),
            _Option(["-offsetfile", "offsetfile"],
                    "transform all features by the offsets given in file",
                    checker_function=check_is_str,
                    filename=True,
                    equate=False,
                    ),
            _Option(["-setsource", "setsource"],
                    "set the source value (2nd column) of each feature",
                    checker_function=check_is_str,
                    equate=False,
                    ),
            _Option(["-typecheck", "typecheck"],
                    ("use an ontology given in an OBO file to validate "
                     "parent-child relationships.  If no argument is given, "
                     "the sofa.obo file from the gtdata/obo_files directory "
                     "is used.  If an argument is given, it is used as an "
                     "OBO filename.  In the case that such a file does not "
                     "exist .obo is added to the argument and loading the "
                     "resulting filename from the gtdata/obo_files "
                     "directory is attempted."),
                    checker_function=check_is_str,
                    filename=True,
                    ),
            _Option(["-xrfcheck", "xrfcheck"],
                    ("check Dbxref and Ontology_term attributes for correct "
                     "syntax according to a abbreviation definition file. "
                     "If no argument is given, the GO.xrf_abbs file from the "
                     "gtdata/xrf_abbr directory is used. If an argument is "
                     "given, it is used as an specific filename for an "
                     "abbreviation file. In the case that such a file does "
                     "not exist, .xrf_abbr is added to the argument and "
                     "loading the resulting filename from the "
                     "gtdata/xrf_abbr directory is attempted."),
                    checker_function=check_is_str,
                    filename=True,
                    ),
            _Switch(["-show", "noshow"],
                    "don't show GFF3 output"
                    ),
            _Switch(["-v", "verbose"],
                    "be verbose"
                    ),
            _Option(["-width", "width"],
                    ("set output width for FASTA sequence printing, "
                     "(0 disables formatting)"),
                    checker_function=check_is_positive_int,
                    equate=False,
                    ),
            _Option(["-o", "outfile"],
                    "redirect output to specified file",
                    checker_function=check_is_str,
                    filename=True,
                    equate=False,
                    ),
            _Switch(["-gzip", "gzip"],
                    "write gzip compressed output file."
                    ),
            _Switch(["-bzip2", "bzip2"],
                    "write bzip2 compressed output file.",
                    ),
            _Switch(["-force", "force"],
                    "force writing to output file",
                    ),
            _Switch(["-help", "help"], "Show help and exit"),
            _Switch(["-version", "version"],
                    "display version information and exit"
                    ),
            _Argument(["infile", "infile2"],
                      "The GFF3 file to operate on.",
                      checker_function=check_is_str,
                      filename=True,
                      is_required=True)
        ]

        super().__init__(cmd, **kwargs)
        return
