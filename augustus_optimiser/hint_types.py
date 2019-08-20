import enum

from typing import Union
from typing import List

# from augustus_optimiser.distributions import DistributionIF
from augustus_optimiser.distributions import AOList
from augustus_optimiser.distributions import FloatConst
from augustus_optimiser.distributions import DistributionIF


class HintType(enum.Enum):

    START = "start"
    STOP = "stop"
    TSS = "tss"
    TTS = "tts"
    ASS = "ass"
    DSS = "dss"
    EXONPART = "exonpart"
    EXON = "exon"
    INTRONPART = "intronpart"
    INTRON = "intron"
    CDSPART = "CDSpart"
    CDS = "CDS"
    UTRPART = "UTRpart"
    UTR = "UTR"
    IRPART = "irpart"
    NONEXONPART = "nonexonpart"
    GENICPART = "genicpart"

    @classmethod
    def parse(cls, string: str) -> "HintType":
        try:
            return STR_TO_HINT[string]
        except KeyError:
            raise ValueError(f"{string} is not a member of {cls.__name__}")

    def __str__(self) -> str:
        return self.value


STR_TO_HINT = {
    "start": HintType.START,
    "stop": HintType.STOP,
    "tss": HintType.TSS,
    "tts": HintType.TTS,
    "ass": HintType.ASS,
    "dss": HintType.DSS,
    "exonpart": HintType.EXONPART,
    "exon": HintType.EXON,
    "intronpart": HintType.INTRONPART,
    "intron": HintType.INTRON,
    "CDSpart": HintType.CDSPART,
    "CDS": HintType.CDS,
    "UTRpart": HintType.UTRPART,
    "UTR": HintType.UTR,
    "irpart": HintType.IRPART,
    "nonexonpart": HintType.NONEXONPART,
    "genicpart": HintType.GENICPART,
}


class HintCell(object):

    def __init__(
        self,
        hint: str,
        boundaries: AOList[Union[float, int]] = AOList([]),
        boni: AOList[Union[float, int]] = AOList([FloatConst(1)]),
    ) -> None:
        """ Creates a weight cell for a specific hint type. """

        self.hint = hint
        self.boundaries = boundaries
        self.boni = boni
        return

    def get(self) -> str:
        """ Get the string representation of a hint cell.

        Examples:
        >>> from augustus_optimiser.distributions import Seq
        >>> from augustus_optimiser.distributions import Range
        >>> from augustus_optimiser.distributions import Uniform
        >>> from augustus_optimiser.distributions import Sample

        >>> h = HintCell('M')
        >>> h.get()
        'M 1 1'

        >>> h = HintCell('M', Seq(5, Range(1, 10)), Seq(6, Uniform(0.5, 1)))
        >>> h.get()
        'M 6 1 4 6 8 10 0.5 0.6 0.7 0.8 0.9 1.0'
        """

        boundaries = self.boundaries.get()
        boni = self.boni.get()
        n_parts = min([len(boundaries) + 1, len(boni)])

        joined_boundaries = " ".join(map(str, boundaries[:n_parts - 1]))
        joined_boni = " ".join(map(str, boni[:n_parts]))

        if n_parts == 1:
            return f"{self.hint} {n_parts} {joined_boni}"
        else:
            return f"{self.hint} {n_parts} {joined_boundaries} {joined_boni}"


class HintRow(object):

    def __init__(
        self,
        kind: HintType,
        bonus: DistributionIF,
        malus: DistributionIF,
        local_malus: DistributionIF,
        hints: List[HintCell],
    ) -> None:
        self.kind = kind
        self.bonus = bonus
        self.malus = malus
        self.local_malus = local_malus
        self.hints = hints
        return

    def get(self) -> str:
        """ Get a string representation of the hint row.

        Examples:
        >>> hr = HintRow(
        ...     HintType.TSS,
        ...     FloatConst(1.0),
        ...     FloatConst(0.5),
        ...     FloatConst(0.1),
        ...     [HintCell('M')],
        ... )
        >>> hr.get()
        'tss 1.0 0.5 0.1 M 1 1'
        """

        bonus = self.bonus.get()
        malus = self.malus.get()
        local_malus = self.local_malus.get()
        hints = " ".join([h.get() for h in self.hints])
        return f"{self.kind} {bonus} {malus} {local_malus} {hints}"
