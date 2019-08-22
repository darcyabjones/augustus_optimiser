"""
"""

import enum

from typing import Union
from typing import Any
from typing import List, Dict
from typing import Sequence, Mapping
from typing import NamedTuple

# from augustus_optimiser.distributions import DistributionIF
from augustus_optimiser.distributions import AOList
from augustus_optimiser.distributions import FloatConst, StrConst
from augustus_optimiser.distributions import Distribution, DistributionIF
from augustus_optimiser.distributions import Choose

from augustus_optimiser.errors import DistValueError


class HintKind(enum.Enum):

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
    def parse(cls, string: str) -> "HintKind":
        try:
            return STR_TO_HINT[string]
        except KeyError:
            raise ValueError(f"{string} is not a member of {cls.__name__}")

    def __str__(self) -> str:
        return self.value

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}.{self.name}"


STR_TO_HINT = {
    "start": HintKind.START,
    "stop": HintKind.STOP,
    "tss": HintKind.TSS,
    "tts": HintKind.TTS,
    "ass": HintKind.ASS,
    "dss": HintKind.DSS,
    "exonpart": HintKind.EXONPART,
    "exon": HintKind.EXON,
    "intronpart": HintKind.INTRONPART,
    "intron": HintKind.INTRON,
    "CDSpart": HintKind.CDSPART,
    "CDS": HintKind.CDS,
    "UTRpart": HintKind.UTRPART,
    "UTR": HintKind.UTR,
    "irpart": HintKind.IRPART,
    "nonexonpart": HintKind.NONEXONPART,
    "genicpart": HintKind.GENICPART,
}


class HintCell(NamedTuple):

    source: str
    boundaries: List[float]
    boni: List[float]

    @classmethod
    def from_obj(cls, d: Mapping[Any, Any]) -> "HintCell":
        # We're going to pretend that this is safe for now.
        # I don't feel like dealing with Optional
        source = d["source"]
        boundaries = d["boundaries"]
        boni = d["boni"]
        return cls(source, boundaries, boni)

    def to_obj(self) -> Dict[str, Union[str, List[float]]]:
        return {
            "source": self.source,
            "boundaries": self.boundaries,
            "boni": self.boni,
        }

    def __str__(self) -> str:
        n_parts = min([len(self.boundaries) + 1, len(self.boni)])

        joined_boundaries = " ".join(map(str, self.boundaries[:n_parts - 1]))
        joined_boni = " ".join(map(str, self.boni[:n_parts]))

        if n_parts == 1:
            return f"{self.source} {n_parts} {joined_boni}"
        else:
            return f"{self.source} {n_parts} {joined_boundaries} {joined_boni}"


class HintCellFactory(object):

    def __init__(
        self,
        source: str,
        boundaries: AOList[Union[float, int]] = AOList([]),
        boni: AOList[Union[float, int]] = AOList([FloatConst(1.0)]),
    ) -> None:
        """ Creates a weight cell for a specific hint type.

        Examples:

        >>> HintCellFactory('M')
        HintCellFactory('M', AOList([]), AOList([FloatConst(1.0)]))
        """

        self.source = source
        self.boundaries = boundaries
        self.boni = boni
        return

    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}("
                f"{repr(self.source)}, "
                f"{repr(self.boundaries)}, "
                f"{repr(self.boni)}"
                ")")

    def get(self) -> HintCell:
        """ Get the string representation of a hint cell.

        Examples:
        >>> from augustus_optimiser.distributions import Seq
        >>> from augustus_optimiser.distributions import Range
        >>> from augustus_optimiser.distributions import Uniform
        >>> from augustus_optimiser.distributions import Sample

        >>> HintCellFactory('M').get()
        HintCell(source='M', boundaries=[], boni=[1.0])

        >>> HintCellFactory(
        ...     'M',
        ...     Seq(5, Range(1, 10)),
        ...     Seq(6, Uniform(0.5, 1.0))
        ... ).get()
        HintCell(source='M',
                 boundaries=[1, 4, 6, 8, 10],
                 boni=[0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        """

        boundaries = self.boundaries.get()
        boni = self.boni.get()

        return HintCell(self.source, boundaries, boni)


HintRowObjType = Dict[
    str,
    Union[
        str,
        float,
        List[Dict[str, Union[str, List[float]]]]
    ]
]


class HintRow(NamedTuple):

    kind: HintKind
    bonus: float
    malus: float
    local_malus: float
    cells: List[HintCell]

    def __str__(self) -> str:
        hints = " ".join(map(str, self.cells))
        return (f"{self.kind} {self.bonus} "
                f"{self.malus} {self.local_malus} {hints}")

    @classmethod
    def from_obj(cls, d: Mapping[Any, Any]) -> "HintRow":
        kind = HintKind.parse(d["kind"])
        cells = [HintCell.from_obj(o) for o in d["cells"]]
        return cls(kind, d["bonus"], d["malus"], d["local_malus"], cells)

    def to_obj(self) -> HintRowObjType:
        return {
            "kind": str(self.kind),
            "bonus": self.bonus,
            "malus": self.malus,
            "local_malus": self.local_malus,
            "cells": [o.to_obj() for o in self.cells]
        }


class HintRowFactory(object):

    def __init__(
        self,
        kind: HintKind,
        bonus: DistributionIF = FloatConst(1.0),
        malus: DistributionIF = FloatConst(1.0),
        local_malus: DistributionIF = FloatConst(1.0),
        cells: Sequence[HintCellFactory] = [
            HintCellFactory("M", boni=AOList([FloatConst(1e100)])),
        ],
    ) -> None:
        """ Represents all hints for a specific type.

        Examples:

        >>> from augustus_optimiser.distributions import Uniform
        >>> HintRowFactory(HintKind.INTRON, bonus=Uniform(1, 2))
        HintRowFactory(HintKind.INTRON,
                Uniform(FloatConst(1.0), FloatConst(2.0)),
                FloatConst(1.0),
                FloatConst(1.0),
                [HintCellFactory('M',
                                 AOList([]),
                                 AOList([FloatConst(1e+100)]))])
        """

        self.kind = kind
        self.bonus = bonus
        self.malus = malus
        self.local_malus = local_malus
        self.cells = list(cells)
        return

    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}("
                f"{repr(self.kind)}, "
                f"{repr(self.bonus)}, "
                f"{repr(self.malus)}, "
                f"{repr(self.local_malus)}, "
                f"{repr(self.cells)}"
                ")")

    def fill_missing(self, sources: List[str]) -> None:
        self.cells = self._add_missing_hints(sources, self.cells)
        return

    @staticmethod
    def _add_missing_hints(
        hints: List[str],
        cells: List[HintCellFactory]
    )-> List[HintCellFactory]:
        """ Adds a default cell for each missing hint type """

        out: List[HintCellFactory] = list()
        cell_dict = {c.source: c for c in cells}

        for hint in hints:
            if hint in cell_dict:
                out.append(cell_dict[hint])
            else:
                out.append(HintCellFactory(hint))
        return out

    def get(self) -> HintRow:
        """ Get a string representation of the hint row.

        Examples:
        >>> HintRowFactory(
        ...     HintKind.TSS,
        ...     FloatConst(1.0),
        ...     FloatConst(0.5),
        ...     FloatConst(0.1),
        ...     [HintCellFactory('M')],
        ... ).get()
        HintRow(kind=HintKind.TSS, bonus=1.0, malus=0.5, local_malus=0.1,
                cells=[HintCell(source='M', boundaries=[], boni=[1.0])])
        """

        bonus = self.bonus.get()
        malus = self.malus.get()
        local_malus = self.local_malus.get()
        hints = [h.get() for h in self.cells]
        return HintRow(self.kind, bonus, malus, local_malus, hints)


class SourceParameter(NamedTuple):

    source: str
    param: str

    def __str__(self) -> str:
        return f"{self.source} {self.param}"

    @classmethod
    def from_obj(cls, d: Mapping[Any, Any]) -> "SourceParameter":
        return cls(d["source"], str(d["param"]))

    def to_obj(self) -> Dict[str, str]:
        return {"source": self.source, "param": self.param}


class SourceParameterFactory(object):

    def __init__(
        self,
        source: str,
        param: Distribution[str] = StrConst("1group1gene")
    ):
        """ Represents a row in the [SOURCE-PARAMETERS] section. """

        self.source = source
        self.param = param

        self.check_valid_dist()
        return

    def __repr__(self):
        return (f"{self.__class__.__name__}"
                f"({repr(self.source)}, {repr(self.param)})")

    def check_valid_dist(self) -> None:
        if isinstance(self.param, StrConst):
            self.check_valid(self.param.value)
        elif isinstance(self.param, Choose):
            for choice in self.param.choices:
                self.check_valid(choice.value)
        return

    @staticmethod
    def check_valid(string: str) -> None:
        valid_strings = ("1group1gene", "individual_liability")
        if string not in valid_strings:
            raise DistValueError(
                f"Received invalid option for source parameter: {string}. "
                f"Source parameters can only be one of {valid_strings}. "
            )
        return

    def get(self) -> SourceParameter:
        """ Gets a string representation of the source parameter row.

        Essentially there are two use cases for this, to specify a constant
        or to randomly select "1group1gene" or "individual_liability".

        Examples:
        >>> SourceParameterFactory("M", StrConst("1group1gene")).get()
        SourceParameter(source='M', param='1group1gene')
        """

        param = self.param.get()
        self.check_valid(param)

        return SourceParameter(self.source, param)


class HintConfig(NamedTuple):

    sources: List[str]
    source_parameters: List[SourceParameter]
    hints: List[HintRow]

    def __str__(self) -> str:
        """ Returns a config file from the object tree.

        Examples:

        >>> hc = HintConfigFactory(["M", "E", "P"])
        >>> print(hc.get())  # doctest: +ELLIPSIS
        # ...
        <BLANKLINE>
        [SOURCES]
        M E P
        <BLANKLINE>
        [SOURCE-PARAMETERS]
        M 1group1gene
        E 1group1gene
        P 1group1gene
        <BLANKLINE>
        [GENERAL]
        start 1.0 1.0 1.0 M 1 1e+100 E 1 1.0 P 1 1.0
        stop 1.0 1.0 1.0 M 1 1e+100 E 1 1.0 P 1 1.0
        tss 1.0 1.0 1.0 M 1 1e+100 E 1 1.0 P 1 1.0
        tts 1.0 1.0 1.0 M 1 1e+100 E 1 1.0 P 1 1.0
        ass 1.0 1.0 1.0 M 1 1e+100 E 1 1.0 P 1 1.0
        dss 1.0 1.0 1.0 M 1 1e+100 E 1 1.0 P 1 1.0
        exonpart 1.0 1.0 1.0 M 1 1e+100 E 1 1.0 P 1 1.0
        exon 1.0 1.0 1.0 M 1 1e+100 E 1 1.0 P 1 1.0
        intronpart 1.0 1.0 1.0 M 1 1e+100 E 1 1.0 P 1 1.0
        intron 1.0 1.0 1.0 M 1 1e+100 E 1 1.0 P 1 1.0
        CDSpart 1.0 1.0 1.0 M 1 1e+100 E 1 1.0 P 1 1.0
        CDS 1.0 1.0 1.0 M 1 1e+100 E 1 1.0 P 1 1.0
        UTRpart 1.0 1.0 1.0 M 1 1e+100 E 1 1.0 P 1 1.0
        UTR 1.0 1.0 1.0 M 1 1e+100 E 1 1.0 P 1 1.0
        irpart 1.0 1.0 1.0 M 1 1e+100 E 1 1.0 P 1 1.0
        nonexonpart 1.0 1.0 1.0 M 1 1e+100 E 1 1.0 P 1 1.0
        genicpart 1.0 1.0 1.0 M 1 1e+100 E 1 1.0 P 1 1.0
        """

        sources_joined = " ".join(self.sources)
        sources = f"[SOURCES]\n{sources_joined}"

        sp_joined = "\n".join(map(str, self.source_parameters))
        source_parameters = f"[SOURCE-PARAMETERS]\n{sp_joined}"

        hints_joined = "\n".join(map(str, self.hints))
        general = f"[GENERAL]\n{hints_joined}"

        return "\n".join([
            "# Augustus extrinsic hints config file.",
            ("# Automatically generated by [augustus_optimiser]"
             "(https://github.com/darcyabjones/augustus_optimiser)"),
            "",
            sources,
            "",
            source_parameters,
            "",
            general,
        ])

    @classmethod
    def from_obj(cls, d: Mapping[Any, Any]) -> "HintConfig":
        """ Construct from a dictionary.

        Examples:
        >>> hc = HintConfigFactory(["M"]).get()
        >>> assert hc == HintConfig.from_obj(hc.to_obj())
        """

        return cls(
            d["sources"],
            [SourceParameter.from_obj(o) for o in d["source_parameters"]],
            [HintRow.from_obj(o) for o in d["hints"]]
        )

    def to_obj(self) -> Dict[str, Union[List[str],
                                        List[Dict[str, str]],
                                        List[HintRowObjType]]]:
        """ Return a dictionary that can be serialised with JSON.

        Examples:

        >>> hc = HintConfigFactory(["M"]).get().to_obj()
        >>> hc["sources"]
        ['M']

        >>> hc["source_parameters"]
        [{'source': 'M', 'param': '1group1gene'}]

        >>> hc["hints"][0]
        {'kind': 'start',
         'bonus': 1.0,
         'malus': 1.0,
         'local_malus': 1.0,
         'cells': [{'source': 'M', 'boundaries': [], 'boni': [1e+100]}]}
        """

        return {
            "sources": self.sources,
            "source_parameters": [
                sp.to_obj()
                for sp
                in self.source_parameters
            ],
            "hints": [h.to_obj() for h in self.hints],
        }


class HintConfigFactory(object):

    def __init__(
        self,
        sources: List[str],
        source_parameters: Sequence[SourceParameterFactory] = [],
        hints: Sequence[HintRowFactory] = [],
    ) -> None:
        """ An object to generate Extrinsic hints config files.

        Parameters:
        sources - A list of the hint sources to include e.g. 'M', 'E' etc.
        source_parameters - A list of SourceParameterFactory objects.
            This tells augustus whether it should use part of a hint group,
            or disregard it entirely if there is a conflict. Any missing
            sources are filled with the default SourceParameterFactory.
        hints - A list of HintRowFactory objects. Missing HintKinds are filled
            with a default row.


        Note that calling `.get()` will modify the `hint_rows` in-place to
        add default values for any missing `hints`. Any hints that are not
        specified in `hints` will be lost.
        """

        self.sources = sources
        self.source_parameters = list(source_parameters)
        self.hints = list(hints)
        return

    def fill_missing(self) -> None:
        """ Adds default values for missing hints in hint rows and
        source parameters.
        """

        self.fill_missing_rows()
        self.fill_missing_source_parameters()
        return

    def fill_missing_source_parameters(self) -> None:
        """ Fill missing hints in source_parameters. """

        sp_dict = {sp.source: sp for sp in self.source_parameters}
        sp_out = list()
        for source in self.sources:
            if source in sp_dict:
                sp_out.append(sp_dict[source])
            else:
                sp_out.append(SourceParameterFactory(source))

        self.source_parameters = sp_out
        return

    def fill_missing_rows(self) -> None:
        """ Fill missing hints and kinds in hint_rows. """

        row_dict = {row.kind: row for row in self.hints}
        out_rows = []

        for kind in HintKind:
            if kind in row_dict:
                row_dict[kind].fill_missing(self.sources)
                out_rows.append(row_dict[kind])
            else:
                new_row = HintRowFactory(kind)
                new_row.fill_missing(self.sources)
                out_rows.append(new_row)

        self.hints = out_rows
        return

    def get(self) -> HintConfig:
        """ Returns a randomly generated config file.

        Examples:

        A config with default parameters and hint sources M.
        The output is very long :(

        >>> hc = HintConfigFactory(["M"])
        >>> x = hc.get()
        >>> assert isinstance(x, HintConfig)
        """

        self.fill_missing()

        source_parameters = [sp.get() for sp in self.source_parameters]

        hints = [row.get() for row in self.hints]
        return HintConfig(self.sources, source_parameters, hints)
