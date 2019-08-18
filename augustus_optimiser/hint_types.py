import enum
from typing import List
from typing import Union
from augustus_optimiser.distributions import PDF, PMF, Categorical


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
        nsteps: Union[int, PMF, Categorical] = 1,
        boundaries: Union[List[float], PDF, Categorical] = [],
        boni: Union[List[float], PDF, PMF, Categorical] = [1]
    ):
        return


class HintRow(object):

    def __init__(
        self,
        kind: HintType,
        bonus: float,
        malus: float,
        local_malus: float,
        custom: List[HintCell],
    ) -> None:
        return
