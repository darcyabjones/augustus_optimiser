import argparse
from copy import deepcopy

from typing import Sequence, Iterator, List, Set, Dict, Optional, Tuple

from intervaltree import Interval, IntervalTree

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation

from gffpal.gff import GFF, GFFRecord

from ao.higher import fmap


EXTRACT_DESC = ("Extract a sub-fasta gff and hints gff given a list of "
                "'good' gene ids to train from.")


def cli_extract(parser):

    parser.add_argument(
        "fasta",
        type=argparse.FileType("r"),
        help="The input fasta to extract training data from.",
    )

    parser.add_argument(
        "gff",
        type=argparse.FileType("r"),
        help="The input gff3 file to extract training data from.",
    )

    parser.add_argument(
        "-i", "--good",
        type=argparse.FileType("r"),
        default=None,
        help=("A file containing a list of gene IDs to fetch, one per line. "
              "If not provided, all genes are used.")
    )

    parser.add_argument(
        "--hints",
        type=argparse.FileType("r"),
        default=None,
        help="Hints for the genome to be coordinate converted."
    )

    parser.add_argument(
        "-g", "--outgff",
        type=argparse.FileType("w"),
        required=True,
        help="Where to write the subset gff3 data out."
    )

    parser.add_argument(
        "-f", "--outfasta",
        type=argparse.FileType("w"),
        required=True,
        help="Where to write the subset fasta data to."
    )

    parser.add_argument(
        "--outhints",
        type=argparse.FileType("w"),
        default=None,
        help="Where to write the subset hints gff to."
    )

    parser.add_argument(
        "-p", "--pad",
        type=int,
        default=100,
        help="Pad the extracted subsections around the genes by this much."
    )

    parser.add_argument(
        "-m", "--merge",
        action="store_true",
        default=False,
        help=(
            "If multiple gene regions to be extracted overlap then merge them "
            "into a single region."
        )
    )

    return


def get_good_ids(handle: Iterator[str]) -> Set[str]:

    out = set()
    for line in handle:
        sline = line.strip()
        out.add(sline)

    return out


def groupby_seqid(gff: GFF) -> Dict[str, List[GFFRecord]]:
    """ Groupby the seqid column of the gff records. """

    from collections import defaultdict
    out: Dict[str, List[GFFRecord]] = defaultdict(list)

    for record in gff:
        out[record.seqid].append(record)

    return out


def gffrecords_to_intervals(features: Sequence[GFFRecord]) -> List[Interval]:
    """ Construct an interval tree from gff records. """

    out = []
    for i, f in enumerate(features):
        if f.start == f.end:
            iv = Interval(f.start, f.start + 1, i)
        elif f.end < f.start:
            iv = Interval(f.end, f.start, i)
        else:
            iv = Interval(f.start, f.end, i)

        out.append(iv)
    return out


def pad_intervals(intervals: Sequence[Interval], pad: int) -> List[Interval]:
    return [Interval(i.begin - pad, i.end + pad, i.data) for i in intervals]


def pad_features(
    features: List[GFFRecord],
    pad: int
) -> Iterator[Tuple[int, int]]:
    for f in features:
        yield f.start - pad, f.end + pad


def subset_features_by_id(
    features: Sequence[GFFRecord],
    target_ids: Set[str],
) -> List[GFFRecord]:
    return [
        f for f in features
        if (f.attributes.id is not None
            and f.attributes.id in target_ids)
    ]


def merge_overlapping(
    features: Sequence[GFFRecord],
    pad: int
) -> Iterator[Tuple[int, int]]:
    """ """

    itree = IntervalTree(
        pad_intervals(
            gffrecords_to_intervals(features),
            pad
        )
    )

    itree.merge_overlaps(strict=False)

    seen: Set[Tuple[int, int]] = set()

    for interval in itree:

        if (interval.begin, interval.end) in seen:
            continue
        else:
            seen.add((interval.begin, interval.end))

        yield interval.begin, interval.end

    return


def shift_gff(gff: GFF, seqid: str, start: int) -> GFF:
    out = list()

    for feature in gff.traverse_children(sort=True):
        f = deepcopy(feature)
        f.start -= start
        f.end -= start
        f.seqid = seqid
        if f.attributes is not None:
            if f.attributes.id is not None:
                f.attributes.id = seqid + "_" + f.attributes.id

            new_parents = []
            for parent in f.attributes.parent:
                new_parent = seqid + "_" + parent
                new_parents.append(new_parent)

            f.attributes.parent = new_parents

        out.append(f)

    return GFF(out)


def make_new_subsequence(
    seqid: str,
    start: int,
    end: int,
    gene_itree: IntervalTree,
    hint_itree: Optional[IntervalTree],
    gene_features: Sequence[GFFRecord],
    hint_features: Optional[Sequence[GFFRecord]],
    seq: SeqRecord
) -> Tuple[str, SeqRecord, GFF, Optional[GFF]]:

    gene_intervals = gene_itree[start: end]

    min_gene_start = min(f.begin for f in gene_intervals) - 10
    max_gene_end = max(f.end for f in gene_intervals) + 10

    start = min([start, min_gene_start])
    if start < 0:
        start = 0

    end = max([end, max_gene_end])
    if end > len(seq):
        end = len(seq)

    if hint_itree is None:
        hint_intervals = None
    else:
        hint_intervals = [
            i for i in hint_itree[start: end]
            if i.begin >= start and i.end <= end
        ]

    name = f"{seqid}:{start}-{end}"

    subseq = FeatureLocation(start, end, 1).extract(seq)
    subseq.id = name
    subseq.name = name
    subseq.description = name

    subgenes = GFF([gene_features[i.data] for i in gene_intervals])
    subgenes_shifted = shift_gff(subgenes, name, start)

    if hint_intervals is None or hint_features is None:
        subhints_shifted = None
    else:
        subhints = GFF([hint_features[i.data] for i in hint_intervals])
        subhints_shifted = shift_gff(subhints, name, start)

    return name, subseq, subgenes_shifted, subhints_shifted


def get_blocks(
    seqs: Dict[str, SeqRecord],
    genes: GFF,
    hints: Optional[GFF],
    pad: int,
    merge: bool,
    target_ids: Optional[Set[str]]
) -> Iterator[Tuple[str, SeqRecord, GFF, Optional[GFF]]]:

    hints_by_seqid = fmap(groupby_seqid, hints)

    for seqid, gene_features in groupby_seqid(
        genes.select_type("gene")
    ).items():
        if seqid not in seqs:
            continue

        gene_itree = IntervalTree(gffrecords_to_intervals(gene_features))
        try:
            hint_itree = fmap(
                IntervalTree,
                fmap(
                    gffrecords_to_intervals,
                    fmap(lambda h: h.get(seqid, []), hints_by_seqid)
                )
            )
        except ValueError as e:
            if hints_by_seqid is not None:
                print(seqid, hints_by_seqid.get(seqid, []))
            raise e

        if target_ids is None:
            target_features: List[GFFRecord] = gene_features
        else:
            target_features = subset_features_by_id(gene_features, target_ids)

        if merge:
            block_iter = merge_overlapping(target_features, pad)
        else:
            block_iter = pad_features(target_features, pad)

        seen: Set[Tuple[int, int]] = set()
        for start, end in block_iter:
            if (start, end) in seen:
                continue
            else:
                seen.add((start, end))

            yield make_new_subsequence(
                seqid,
                start,
                end,
                gene_itree,
                hint_itree,
                gene_features,
                fmap(lambda h: h.get(seqid, []), hints_by_seqid),
                seqs[seqid],
            )

    return


def run_extract(args) -> None:

    if args.good is not None:
        good_ids: Optional[Set[str]] = get_good_ids(args.good)
    else:
        good_ids = None

    in_gff = GFF.parse(args.gff)
    in_fasta = SeqIO.to_dict(SeqIO.parse(args.fasta, format="fasta"))

    if args.hints is not None:
        in_hints = GFF.parse(args.hints)
    else:
        in_hints = None

    block_iter = get_blocks(
        in_fasta,
        in_gff,
        in_hints,
        args.pad,
        args.merge,
        good_ids
    )

    print("##gff-version 3", file=args.outgff)

    if args.outhints is not None:
        print("##gff-version 3", file=args.outhints)

    for name, seq, genes, hints in block_iter:
        SeqIO.write(seq, args.outfasta, format="fasta")

        sr = f"##sequence-region   {name} 1 {len(seq)}"
        print(sr, file=args.outgff)
        print(genes, file=args.outgff)

        if (args.outhints is not None
                and hints is not None
                and len(hints.inner) > 0):
            print(sr, file=args.outhints)
            print(hints, file=args.outhints)

    return
