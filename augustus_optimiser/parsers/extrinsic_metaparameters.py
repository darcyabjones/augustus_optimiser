from typing import Any
from typing import Dict, List, Tuple
from typing import Union

from augustus_optimiser.distributions import Distribution, DistributionIF
from augustus_optimiser.distributions import FloatConst
from augustus_optimiser.distributions import AOList

from augustus_optimiser.hints import HintConfigFactory
from augustus_optimiser.hints import SourceParameterFactory
from augustus_optimiser.hints import HintRowFactory, HintCellFactory
from augustus_optimiser.hints import HintKind

from augustus_optimiser.errors import DistParseError, ConfigParseError
from augustus_optimiser.parsers.distribution_functions import (
    parse_float_dist, parse_float_list, parse_str_dist
)


def parser(config: Dict[Any, Any]) -> HintConfigFactory:
    """ Parses a dictionary representation of the config file. """

    sources = parse_sources(config.get("sources", None))
    source_parameters = parse_source_parameters(
        config.get("source_parameters", None)
    )

    weights: Dict[HintKind, Dict[str, DistributionIF]] = dict()
    for kind, column, dist in parse_weights(config.get("weights", None)):
        if kind in weights:
            weights[kind][column] = dist
        else:
            weights[kind] = {column: dist}

    source_weights = parse_source_weights(config.get("source_weights", None))

    hint_rows = list()
    for kind in HintKind:
        kind_weights = weights.get(kind, {})
        kind_source_weights = [v for k, v in source_weights if k == kind]
        this_row = HintRowFactory(
            kind,
            bonus=kind_weights.get("bonus", FloatConst(1)),
            malus=kind_weights.get("malus", FloatConst(1)),
            local_malus=kind_weights.get("local_malus", FloatConst(1)),
            cells=kind_source_weights
        )
        hint_rows.append(this_row)

    return HintConfigFactory(sources, source_parameters, hint_rows)


def get_from_object(key: str, obj: Dict, table: str) -> Any:
    """ Helper to get a field from an object or raise an error. """

    val = obj.get(key, None)
    if val is None:
        raise ConfigParseError(
            f"The required field '{key}' was missing from object {obj} in the "
            "'{table}'."
        )

    return val


def get_kind_from_object(obj: Dict, table: Any) -> HintKind:
    kind_str = get_from_object("kind", obj, table)

    if not isinstance(kind_str, str):
        raise ConfigParseError(
            f"Error while parsing '{table}' object {obj}. "
            "The field 'kind' must be a string."
        )

    try:
        kind = HintKind.parse(kind_str)
    except ValueError:
        valid_kinds = str(list(map(str, HintKind)))
        raise ConfigParseError(
            f"Error while parsing '{table}' object {obj}. "
            "The field 'kind' has an invalid value. "
            f"Must be one of {valid_kinds}."
        )

    return kind


def try_catch_parse_dist(
    dist_obj: Any,
    kind: HintKind,
    field: str,
    table: str
) -> DistributionIF:
    try:
        this_dist = parse_float_dist(dist_obj)
    except DistParseError as e:
        raise ConfigParseError(
            f"Error while parsing distribution for '{kind}' '{field}' "
            f"in '{table}'. Input was {dist_obj}. "
            f"{e.msg}"
        )
    return this_dist


def try_catch_parse_list(
    dist_obj: Any,
    kind: HintKind,
    field: str,
    table: str
) -> AOList[Union[int, float]]:
    try:
        this_dist = parse_float_list(dist_obj)
    except DistParseError as e:
        raise ConfigParseError(
            f"Error while parsing distribution for '{kind}' '{field}' "
            f"in '{table}'. Input was {dist_obj}. "
            f"{e.msg}"
        )
    return this_dist


def parse_source_weights(
    weights: Any
) -> List[Tuple[HintKind, HintCellFactory]]:
    """ Parse the weight table from the config file. """

    if weights is None or (isinstance(weights, list) and len(weights) == 0):
        # We don't have any weights to optimize.
        return []

    elif not isinstance(weights, list):
        raise ConfigParseError(
            "The 'source_weights' section of the config file must be a "
            f"list of objects. You provided '{weights}'."
        )

    offensive_elements = [w for w in weights if not isinstance(w, dict)]
    if len(offensive_elements) > 0:
        raise ConfigParseError(
            "The 'source_weights' section of the config file must be a "
            "list of objects. The following elements are not objects: "
            f"{offensive_elements}."
        )

    out = list()
    for weight in weights:
        kind = get_kind_from_object(weight, "weights")
        source = get_from_object("source", weight, "weights")

        if "weight" in weight:
            if ("weights" in weight) or ("boundaries" in weight):
                raise ConfigParseError(
                    f"Error while parsing 'source_weights' for {kind} "
                    f"{source}. "
                    "The 'weight' and the 'weights'/'boundaries' fields are "
                    "mutually exclusive. Please remove the fields that you "
                    "don't intend to use."
                )

            this_dist: Distribution = try_catch_parse_dist(
                weight["weight"],
                kind,
                "weight",
                "source_weights"
            )

            this_list_dist: AOList[Union[int, float]] = AOList([this_dist])
            out.append((kind, HintCellFactory(source, boni=this_list_dist)))

        elif ("weights" in weight) and ("boundaries" in weight):
            weights_list = weight["weights"]
            weights_dists = try_catch_parse_list(
                weights_list,
                kind,
                "weights",
                "source_weights"
            )

            boundaries_list = weight["boundaries"]
            boundaries_dists = try_catch_parse_list(
                boundaries_list,
                kind,
                "boundaries",
                "source_weights"
            )
            out.append((
                kind,
                HintCellFactory(source, boundaries_dists, weights_dists)
            ))

        elif ("weights" in weight) or ("boundaries" in weight):
            raise ConfigParseError(
                f"Error while parsing '{kind}' '{source}' "
                "in 'source_weights'. "
                "When using the 'weights'/'boundaries' fields, both must be "
                "provided. Please add the missing field or use the 'weight' "
                "field if you are providing a single distribution."
            )
        else:
            raise ConfigParseError(
                "Error while parsing 'source_weights'. "
                "Either the 'weight' or both the 'weights' and 'boundaries' "
                "fields must be provided."
            )

    return out


def parse_weights(weights: Any) -> List[Tuple[HintKind, str, Distribution]]:
    """ Parse the weight table from the config file.

    Examples:
    >>> parse_weights([])
    []

    >>> parse_weights(None)
    []

    >>> parse_weights([
    ...     {"kind": "exon", "malus": 0.9, "bonus": ["uniform", 1, 2]},
    ...     {"kind": "exonpart", "malus": 0.99, "local_malus": ["beta", 1, 1]},
    ... ])
    [(HintKind.EXON, 'bonus', Uniform(FloatConst(1.0), FloatConst(2.0))),
    (HintKind.EXON, 'malus', FloatConst(0.9)),
    (HintKind.EXONPART, 'malus', FloatConst(0.99)),
    (HintKind.EXONPART, 'local_malus', Beta(FloatConst(1.0), FloatConst(1.0)))]
    """

    if weights is None or (isinstance(weights, list) and len(weights) == 0):
        # We don't have any weights to optimize.
        return []

    elif not isinstance(weights, list):
        raise ConfigParseError(
            "The 'weights' section of the config file must be a list of "
            f"objects. You provided '{weights}'."
        )

    offensive_elements = [w for w in weights if not isinstance(w, dict)]
    if len(offensive_elements) > 0:
        raise ConfigParseError(
            "The 'weights' section of the config file must be a list of "
            "objects. The following elements are not objects: "
            f"{offensive_elements}."
        )

    out = list()
    for weight in weights:
        kind = get_kind_from_object(weight, "weights")

        for column in ("bonus", "malus", "local_malus"):
            if column not in weight:
                continue

            this_dist = try_catch_parse_dist(
                weight[column],
                kind,
                column,
                "weights"
            )
            out.append((kind, column, this_dist))

    return out


def parse_sources(sources: Any) -> List[str]:
    """ An interface that checks that the right types were given.

    Examples:
    >>> parse_sources(["one", "two", "three"])
    ['one', 'two', 'three']

    >>> parse_sources([])
    Traceback (most recent call last):
        ...
    ConfigParseError: ...

    >>> parse_sources("M")
    Traceback (most recent call last):
        ...
    ConfigParseError: ...

    >>> parse_sources(["one", 2, "three"])
    Traceback (most recent call last):
        ...
    ConfigParseError: ...
    """

    if sources is None or (isinstance(sources, list) and len(sources) == 0):
        raise ConfigParseError(
            "The 'sources' section of the config file was missing or empty. "
            "We can't optimise hints when there are no hints configured."
        )

    elif not isinstance(sources, list):
        raise ConfigParseError(
            "The 'sources' section of the config file must be a list of "
            f"strings. You provided {repr(sources)}."
        )

    offensive_elements = [v for v in sources if not isinstance(v, str)]

    if len(offensive_elements) > 0:
        raise ConfigParseError(
            "The 'sources' section of the config file must be a list "
            "of strings. The following elements are not strings: "
            f"{offensive_elements}. Consider quoting them if they are "
            "the sources you want."
        )

    # We know that they are all strings because offensive_elements is empty.
    return sources


def parse_source_parameters(sps: Any) -> List[SourceParameterFactory]:
    """ Checks that the right types were given and parses them to py objects.

    This section is optional so None or empty lists are ok.

    Expected input is of the form:
    [{"source": str, "parameter": param}]"

    Examples:
    >>> parse_source_parameters([{"source": "M", "parameter": "1group1gene"}])
    [SourceParameterFactory('M', StrConst('1group1gene'))]

    >>> parse_source_parameters([])
    []

    >>> parse_source_parameters(None)
    []

    >>> parse_source_parameters({"source": "M", "parameter": "1group1gene"})
    Traceback (most recent call last):
        ...
    ConfigParseError: ...

    >>> parse_source_parameters([
    ...     {"source": "M",
    ...      "parameter": ["choose", ["1group1gene", "individual_liability"]]}
    ... ])
    [SourceParameterFactory('M', ChooseS([StrConst('1group1gene'),
                                   StrConst('individual_liability')]))]
    >>> parse_source_parameters([
    ...     {"source": "M", "parameter": "1group1gene"},
    ...     {"source": "E", "parameter": "individual_liability"}
    ... ])
    [SourceParameterFactory('M', StrConst('1group1gene')),
     SourceParameterFactory('E', StrConst('individual_liability'))]
    """

    if sps is None:
        return []
    elif not isinstance(sps, list):
        raise ConfigParseError(
            "The 'source_parameters' section of the config file must be a "
            "list of objects. "
            "You provided {sps}."
        )

    offensive_elements = [sp for sp in sps if not isinstance(sp, dict)]
    if len(offensive_elements) > 0:
        raise ConfigParseError(
            "The 'source_parameters' section of the config file must be a "
            "list of objects. "
            f"The following elements are not objects: {offensive_elements}."
        )

    out = list()
    for sp in sps:

        source = get_from_object("source", sp, "source_parameters")
        param = get_from_object("parameter", sp, "source_parameters")

        try:
            param_dist = parse_str_dist(param)
        except DistParseError as e:
            raise ConfigParseError(
                f"Error while parsing 'source_parameter' for {source}. "
                f"{e.msg}"
            )
        out.append(SourceParameterFactory(source, param_dist))

    return out
