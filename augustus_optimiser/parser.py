from typing import Any
from typing import Sequence
from typing import Optional
from typing import Callable
from typing import Union
from typing import List


from augustus_optimiser.distributions import Distribution
from augustus_optimiser.distributions import DistributionIF
from augustus_optimiser.distributions import PDF, PMF, Categorical
from augustus_optimiser.distributions import FloatConst, IntConst, StrConst
from augustus_optimiser.distributions import Uniform, ExpUniform
from augustus_optimiser.distributions import Beta, Chi, Gamma, Normal
from augustus_optimiser.distributions import ChooseI, ChooseF, ChooseS
from augustus_optimiser.distributions import Range
from augustus_optimiser.distributions import Binomial, NegBinomial
from augustus_optimiser.distributions import HyperGeometric
from augustus_optimiser.distributions import Trunc
from augustus_optimiser.distributions import AOList

from augustus_optimiser.errors import DistParseError


def str_to_int_parser(
    fn: str
) -> Optional[Callable[[List[Any]], Distribution[int]]]:
    """ Takes a string function name, returns the parser function """

    str_to_int = {
        "range": parse_range,
        "choose": parse_choose_int,
        "binomial": parse_binomial,
        "nbinomial": parse_neg_binomial,
        "hypergeom": parse_hypergeom,
        "trunc": parse_trunc_int,
    }
    return str_to_int.get(fn, None)


def str_to_float_parser(
    fn: str
) -> Optional[Callable[[List[Any]], DistributionIF]]:
    """ Takes a string function name, returns the parser function """

    str_to_float = {
        "uniform": parse_uniform,
        "expuniform": parse_expuniform,
        "beta": parse_beta,
        "chi": parse_chi,
        "gamma": parse_gamma,
        "normal": parse_normal,
        "choose": parse_choose_float,
        "trunc": parse_trunc_float,
    }

    float_fn = str_to_float.get(fn, None)
    if float_fn is not None:
        return float_fn

    int_fn = str_to_int_parser(fn)
    if int_fn is None:
        return None

    return int_fn


def str_to_str_parser(
    fn: str
) -> Optional[Callable[[List[Any]], Distribution[str]]]:
    """ Takes a string function name, returns the parser function """

    if fn == "choose":
        return parse_choose_str
    else:
        return None


# Distribution parsers

def parse_string_as_int(elem: Union[int, str]) -> Distribution[int]:
    try:
        parsed_elem: Distribution[int] = IntConst(int(elem))
    except ValueError:
        raise DistParseError(f"Could not parse value '{elem}' as int.")

    return parsed_elem


def parse_list_in_int_dist(elem: List[Any]) -> Distribution[int]:
    """ Parse a list while parsing a distribution.

    This could be a function or an actual list.
    """

    if len(elem) == 0:
        raise DistParseError("Encountered an empty list")

    elif isinstance(elem[0], str):
        func = str_to_int_parser(elem[0])
        if func is None:
            raise DistParseError(
                f"Recieved a list or invalid function {elem}."
            )
        else:
            parsed_elem = func(elem[1:])

    else:
        raise DistParseError(
            f"Received a list but cannot take a list here: {elem}"
        )
    return parsed_elem


def parse_int_dist(elem: Any) -> Distribution[int]:
    """ Parse an integer distribution specification.

    Examples:

    >>> parse_float_dist(['range', 1, 2])
    Range(IntConst(1), IntConst(2))
    >>> parse_int_dist(1)
    IntConst(1)
    >>> parse_int_dist(['choose', [3, 5, 7]])
    ChooseI([IntConst(3), IntConst(5), IntConst(7)])
    >>> parse_int_dist([1, 2, 3])
    Traceback (most recent call last):
        ...
    DistParseError: Received a list but cannot take a list here: {...}.
    """

    if isinstance(elem, (int, str)):
        parsed_elem = parse_string_as_int(elem)

    elif isinstance(elem, list):
        parsed_elem = parse_list_in_int_dist(elem)

    elif isinstance(elem, Distribution):
        parsed_elem = elem

    elif isinstance(elem, float):
        raise DistParseError(
            "This distribution doesn't support float values. "
            "Please convert '{elem}' to an integer."
        )

    else:
        raise DistParseError(
            "Recieved invalid input: '{repr(elem)}'."
        )

    if isinstance(parsed_elem, PMF):
        return parsed_elem
    else:
        raise DistParseError(
            f"Could not parse value '{elem}' as an integer distribution."
        )


def parse_string_as_float(elem: Union[int, float, str]) -> DistributionIF:
    try:
        parsed_elem: DistributionIF = FloatConst(float(elem))
    except ValueError:
        raise DistParseError(f"Could not parse value '{elem}' as float.")
    return parsed_elem


def parse_list_in_float_dist(elem: List[Any]) -> DistributionIF:
    if len(elem) == 0:
        raise DistParseError("Encountered an empty list")

    elif isinstance(elem[0], str):
        func = str_to_float_parser(elem[0])
        if func is None:
            raise DistParseError(
                f"Recieved a list or invalid function {elem}."
            )
        else:
            parsed_elem = func(elem[1:])
    else:
        raise DistParseError(
            f"Received a list but cannot take a list here: {elem}"
        )
    return parsed_elem


def parse_float_dist(elem: Any) -> DistributionIF:
    """ Parse a float distribution specification.

    Examples:

    >>> parse_float_dist(['uniform', 1, 2])
    Uniform(FloatConst(1.0), FloatConst(2.0))
    >>> parse_float_dist(['uniform', 1, ['uniform', 3, 4.5]])
    Uniform(FloatConst(1.0), Uniform(FloatConst(3.0), FloatConst(4.5)))
    >>> parse_float_dist(['range', 1, 3])
    Range(IntConst(1), IntConst(3))
    >>> parse_float_dist(['uniform', 1, ['range', 3, 7]])
    Uniform(FloatConst(1.0), Range(IntConst(3), IntConst(7)))
    >>> parse_float_dist(['choose', [3, 5, 7]])
    ChooseF([FloatConst(3.0), FloatConst(5.0), FloatConst(7.0)])
    >>> parse_float_dist([1, 2, 3])
    Traceback (most recent call last):
        ...
    DistParseError: Received a list but cannot take a list here: {...}.
    """

    if isinstance(elem, (int, float, str)):
        parsed_elem = parse_string_as_float(elem)

    elif isinstance(elem, list):
        parsed_elem = parse_list_in_float_dist(elem)

    elif isinstance(elem, Distribution):
        parsed_elem = elem

    else:
        raise DistParseError(
            "Recieved invalid input: '{repr(elem)}'."
        )

    if isinstance(parsed_elem, (PDF, PMF)):
        return parsed_elem
    else:
        raise DistParseError(
            f"Could not parse value '{elem}' as a float distribution."
        )


def parse_list_in_str_dist(elem: List[Any]) -> Distribution[str]:
    if len(elem) == 0:
        raise DistParseError("Encountered an empty list")

    elif isinstance(elem[0], str):
        func = str_to_str_parser(elem[0])
        if func is None:
            raise DistParseError(
                f"Recieved a list or invalid function {elem}."
            )
        else:
            parsed_elem = func(elem[1:])

    else:
        raise DistParseError(
            f"Received a list but cannot take a list here: {elem}"
        )
    return parsed_elem


def parse_str_dist(elem: Any) -> Distribution[str]:
    if isinstance(elem, (int, float, str)):
        parsed_elem: Distribution[str] = StrConst(str(elem))

    elif isinstance(elem, list):
        parsed_elem = parse_list_in_str_dist(elem)

    elif isinstance(elem, Distribution):
        parsed_elem = elem

    else:
        raise DistParseError(
            "Recieved invalid input: '{repr(elem)}'."
        )

    if isinstance(parsed_elem, Categorical):
        return parsed_elem

    else:
        raise DistParseError(
            f"Could not parse value '{elem}' as str distribution."
        )


# List parsers
#
# These are types that return a list of distributions, rather than converting
# a list into a single distribution.

def parse_int_list(elems: Sequence[Any]) -> AOList[int]:
    """ Parse a list of integer distributions. """

    try:
        return AOList([parse_int_dist(e) for e in elems])
    except DistParseError as e:
        raise DistParseError(
            f"Error while parsing integer list {repr(elems)}. {e.msg}"
        )


def parse_float_list(elems: Sequence[Any]) -> AOList[Union[int, float]]:
    """ Parse a list of float distributions.

    Examples:
    >>> parse_float_list([1, 2, 3])
    AOList([FloatConst(1.0), FloatConst(2.0), FloatConst(3.0)])
    >>> parse_float_list([2, ['uniform', 3, 4.6]])
    AOList([FloatConst(2.0), Uniform(FloatConst(3.0), FloatConst(4.6))])
    """

    try:
        # I'm being a bit cheeky with type hints here.
        # Return from parse_float_dist is DistributionIF
        # The type system isn't expressive enough to destructure to
        # AOList[Union[int, float]].
        comp: Sequence[Distribution] = [parse_float_dist(e) for e in elems]
        return AOList(comp)
    except DistParseError as e:
        raise DistParseError(
            f"Error while parsing float list {repr(elems)}. {e.msg}"
        )


def parse_str_list(elems: Sequence[Any]) -> AOList[str]:
    """ Parse a list of string distributions. """

    try:
        return AOList([parse_str_dist(e) for e in elems])
    except DistParseError as e:
        raise DistParseError(
            f"Error while parsing str list {repr(elems)}. {e.msg}"
        )


# Function parsers
# Parse a distribution given a list of parameters.

def check_float_param(
    index: int,
    elems: Sequence[Any],
    name: str
) -> DistributionIF:
    """ Does error handling and conversion of float parameters. """

    try:
        param = parse_float_dist(elems[index])
    except IndexError:
        raise DistParseError(
            f"Missing <{name}> parameter for pdf.",
        )
    except ValueError:
        raise DistParseError(
            f"Could not parse <{name}> parameter '{elems[index]}' as float.",
        )
    return param


def check_int_param(
    index: int,
    elems: Sequence[Any],
    name: str
) -> Distribution[int]:
    """ Does error handling and conversion of int parameters. """

    try:
        param = parse_int_dist(elems[index])
    except IndexError:
        raise DistParseError(
            f"Missing <{name}> parameter for pmf.",
        )
    except ValueError:
        raise DistParseError(
            f"Could not parse <{name}> parameter '{elems[index]}' as integer.",
        )
    return param


def parse_uniform(elems: Sequence[Any]) -> Uniform:
    """ Parse a list definition as a uniform distribution.

    Examples:

    >>> u = parse_uniform([1, 3.5])
    >>> u.min
    FloatConst(1.0)
    >>> u.max
    FloatConst(3.5)
    >>> parse_uniform(["one", 2])
    Traceback (most recent call last):
        ...
    DistParseError: Could not parse min value 'one' as float.
    """

    min_ = check_float_param(0, elems, "min")
    max_ = check_float_param(1, elems, "max")

    # TODO, figure out how to assert this now that either could be any
    # distribution
    # if min_ >= max_:
    #     raise DistParseError(f"The <min> bound must be smaller than <max>")

    return Uniform(min_, max_)


def parse_expuniform(elems: Sequence[Any]) -> ExpUniform:
    """ Parse a list definition as an exponential uniform distribution.

    Examples:

    >>> u = parse_expuniform([10, 1, 3.5])
    >>> u.base
    FloatConst(10.0)
    >>> u.min
    FloatConst(1.0)
    >>> u.max
    FloatConst(3.5)
    >>> parse_expuniform(["one", 2])
    Traceback (most recent call last):
        ...
    DistParseError: Could not parse base value 'one' as float.
    """

    base = check_float_param(0, elems, "base")
    min_ = check_float_param(1, elems, "min")
    max_ = check_float_param(2, elems, "max")

    # TODO, figure out how to assert this now that either could be any
    # distribution
    # if min_ >= max_:
    #     raise DistParseError(f"The <min> bound must be smaller than <max>")

    return ExpUniform(base, min_, max_)


def parse_beta(elems: Sequence[Any]) -> Beta:
    """ Parse a list definition as a beta distribution.

    Examples:

    >>> u = parse_beta([2.31, 0.627])
    >>> u.a
    FloatConst(2.31)

    >>> u.b
    FloatConst(0.627)

    >>> parse_beta(["one", 2])
    Traceback (most recent call last):
        ...
    DistParseError: Could not parse a value 'one' as float.
    """

    a = check_float_param(0, elems, "a")
    b = check_float_param(1, elems, "b")

    # TODO, figure out how to assert this now that either could be any
    # distribution
    # if min_ >= max_:
    #     raise DistParseError(f"The <min> bound must be smaller than <max>")

    return Beta(a, b)


def parse_chi(elems: Sequence[Any]) -> Chi:
    """ Parse a list definition as a chi distribution.

    Examples:

    >>> u = parse_chi([1, 2.31, 0.627])
    >>> u.df
    FloatConst(1.0)

    >>> u.loc
    FloatConst(2.31)

    >>> u.scale
    FloatConst(0.627)

    >>> parse_chi(["one", 2])
    Traceback (most recent call last):
        ...
    DistParseError: Could not parse a value 'one' as float.
    """

    df = check_float_param(0, elems, "df")
    loc = check_float_param(1, elems, "loc")
    scale = check_float_param(2, elems, "scale")

    # TODO, figure out how to assert this now that either could be any
    # distribution
    # if min_ >= max_:
    #     raise DistParseError(f"The <min> bound must be smaller than <max>")

    return Chi(df, loc, scale)


def parse_gamma(elems: Sequence[Any]) -> Gamma:
    """ Parse a list definition as a chi distribution definition.

    Examples:

    >>> u = parse_gamma([1, 2.31, 0.627])
    >>> u.a
    FloatConst(1.0)

    >>> u.loc
    FloatConst(2.31)

    >>> u.scale
    FloatConst(0.627)

    >>> parse_gamma(["one", 2])
    Traceback (most recent call last):
        ...
    DistParseError: Could not parse a value 'one' as float.
    """

    a = check_float_param(0, elems, "a")
    loc = check_float_param(1, elems, "loc")
    scale = check_float_param(2, elems, "scale")

    return Gamma(a, loc, scale)


def parse_normal(elems: Sequence[Any]) -> Normal:
    """ Parse a list definition as a normal distribution definition.

    Examples:

    >>> u = parse_normal([1, 2.31])
    >>> u.loc
    FloatConst(1.0)

    >>> u.scale
    FloatConst(2.31)

    >>> parse_gamma(["one", 2])
    Traceback (most recent call last):
        ...
    DistParseError: Could not parse a value 'one' as float.
    """

    loc = check_float_param(0, elems, "loc")
    scale = check_float_param(1, elems, "scale")

    return Normal(loc, scale)


def parse_range(elems: Sequence[Any]) -> Range:
    """ Parse a list definition as a range distribution definition.

    Examples:

    >>> r = parse_range([1, 3])
    >>> r.min
    IntConst(1)

    >>> r.max
    IntConst(3)

    >>> r = parse_range([1, 2, 3])
    >>> r.max
    IntConst(2)
    """

    min_ = check_int_param(0, elems, "min")
    max_ = check_int_param(1, elems, "max")

    # if min_ >= max_:
    #     raise DistParseError(f"The <min> bound must be smaller than <max>")

    return Range(min_, max_)


def parse_binomial(elems: Sequence[Any]) -> Binomial:
    """ Parse a list definition as a binomial distribution definition.

    Examples:

    >>> r = parse_binomial([5, 0.4, 2])
    >>> r.n
    IntConst(5)

    >>> r.p
    FloatConst(0.4)

    >>> r.loc
    IntConst(2)
    """

    n = check_int_param(0, elems, "n")
    p = check_float_param(1, elems, "p")

    if not isinstance(p, PDF):
        raise DistParseError(
            "Parameter <p> to binomial must be a float between 0 and 1. "
            "You provided '{elems[1]}'."
        )

    loc = check_int_param(2, elems, "loc")

    return Binomial(n, p, loc)


def parse_neg_binomial(elems: Sequence[Any]) -> NegBinomial:
    """ Parse a list definition as a negative binomial distribution definition.

    Examples:

    >>> r = parse_neg_binomial([5, 0.4, 2])
    >>> r.n
    IntConst(5)

    >>> r.p
    FloatConst(0.4)

    >>> r.loc
    IntConst(2)
    """

    n = check_int_param(0, elems, "n")
    p = check_float_param(1, elems, "p")

    if not isinstance(p, PDF):
        raise DistParseError(
            "Parameter <p> to nbinomial must be a float between 0 and 1. "
            "You provided '{elems[1]}'."
        )

    loc = check_int_param(2, elems, "loc")

    return NegBinomial(n, p, loc)


def parse_hypergeom(elems: Sequence[Any]) -> HyperGeometric:
    """ Parse a list definition as a hypergeometric distribution definition.

    Examples:

    >>> r = parse_hypergeom([5, 8, 12, 0])
    >>> r.n_success
    IntConst(5)

    >>> r.n
    IntConst(8)

    >>> r.population
    IntConst(12)

    >>> r.loc
    IntConst(0)
    """

    n_success = check_int_param(0, elems, "n_success")
    n = check_int_param(1, elems, "n")
    population = check_int_param(2, elems, "population")
    loc = check_int_param(3, elems, "loc")

    return HyperGeometric(n_success, n, population, loc)


def parse_trunc_int(elems: Sequence[Any]) -> Trunc[int]:
    """ Parse a list definition as a truncated distribution definition.

    Examples:

    >>> r = parse_trunc_int([5, 10, ["range", 1, 20]])
    >>> r.min
    5

    >>> r.max
    10

    >>> r.dist
    Range(IntConst(1), IntConst(20))
    """

    min_ = check_int_param(0, elems, "min")
    if not isinstance(min_, IntConst):
        raise DistParseError(
            "The <min> parameter to trunc must be a single integer. "
            "It cannot be a distribution. "
            f"You provided {elems[0]}."
        )
    max_ = check_int_param(1, elems, "max")
    if not isinstance(max_, IntConst):
        raise DistParseError(
            "The <max> parameter to trunc must be a single integer. "
            "It cannot be a distribution. "
            f"You provided {elems[1]}."
        )
    dist = check_int_param(2, elems, "dist")

    return Trunc(min_.value, max_.value, dist)


def parse_trunc_float(elems: Sequence[Any]) -> Trunc[Union[float, int]]:
    """ Parse a list definition as a truncated distribution definition.

    Examples:

    >>> r = parse_trunc_float([5, 10, ["uniform", 1, 20]])
    >>> r.min
    5.0

    >>> r.max
    10.0

    >>> r.dist
    Uniform(FloatConst(1.0), FloatConst(20.0))
    """

    min_ = check_float_param(0, elems, "min")
    if not isinstance(min_, FloatConst):
        raise DistParseError(
            "The <min> parameter to trunc must be a single float. "
            "It cannot be a distribution. "
            f"You provided {elems[0]}."
        )

    max_ = check_float_param(1, elems, "max")
    if not isinstance(max_, FloatConst):
        raise DistParseError(
            "The <max> parameter to trunc must be a single float. "
            "It cannot be a distribution. "
            f"You provided {elems[1]}."
        )

    # I'm being a bit cheeky with type hints here.
    # Return from check float param is DistributionIF
    # The type system isn't expressive enough to destructure to
    # Trunc[Union[int, float]].
    dist: Distribution[Any] = check_float_param(2, elems, "dist")
    return Trunc(min_.value, max_.value, dist)


def parse_choose_int(elems: Sequence[Any]) -> ChooseI:
    """ Parse choose from a list representation.

    Examples:
    >>> c = parse_choose_int([[1, 2, ['range', 5, 10]]])
    >>> c.choices
    [IntConst(1), IntConst(2), Range(IntConst(5), IntConst(10))]
    """

    if len(elems) == 0:
        raise DistParseError("No options specified for choose to use.")
    elif len(elems) > 1:
        raise DistParseError(f"Choose takes only one argument: {elems}")
    else:
        parsed_elems = parse_int_list(elems[0])

    return ChooseI(parsed_elems.elems)


def parse_choose_float(elems: Sequence[Any]) -> ChooseF:
    """ Parse choose from a list representation.

    Examples:
    >>> c = parse_choose_float([[2, ['uniform', 5, 10]]])
    >>> c.choices
    [FloatConst(2.0), Uniform(FloatConst(5.0), FloatConst(10.0))]
    """

    if len(elems) == 0:
        raise DistParseError("No options specified for choose to use.")
    elif len(elems) > 1:
        raise DistParseError(f"Choose takes only one argument: {elems}")
    else:
        parsed_elems = parse_float_list(elems[0])

    return ChooseF(parsed_elems.elems)


def parse_choose_str(elems: Sequence[Any]) -> ChooseS:
    """ Parse choose from a list representation.

    Examples:
    >>> c = parse_choose_str([[2, "two"]])
    >>> c.choices
    [StrConst('2'), StrConst('two')]
    """

    if len(elems) == 0:
        raise DistParseError("No options specified for choose to use.")
    elif len(elems) > 1:
        raise DistParseError(f"Choose takes only one argument: {elems}")
    else:
        parsed_elems = parse_str_list(elems[0])

    return ChooseS(parsed_elems.elems)
