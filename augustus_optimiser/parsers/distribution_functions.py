from typing import Any
from typing import Sequence
from typing import Optional
from typing import Callable
from typing import Union
from typing import List
from typing import TypeVar

from augustus_optimiser.distributions import Distribution
from augustus_optimiser.distributions import DistributionIF
from augustus_optimiser.distributions import AOList
from augustus_optimiser import distributions as dist

from augustus_optimiser.errors import DistParseError


T = TypeVar('T')


def int_dist_fn_finder(
    fn: str
) -> Optional[Callable[[Sequence[Any]], Distribution[int]]]:
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


def float_dist_fn_finder(
    fn: str
) -> Optional[Callable[[Sequence[Any]], DistributionIF]]:
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

    int_fn = int_dist_fn_finder(fn)
    if int_fn is None:
        return None

    return int_fn


def str_dist_fn_finder(
    fn: str
) -> Optional[Callable[[Sequence[Any]], Distribution[str]]]:
    """ Takes a string function name, returns the parser function """

    if fn == "choose":
        return parse_choose_str
    else:
        return None


def float_list_fn_finder(
    fn: str
) -> Optional[Callable[[Sequence[Any]], AOList[Union[float, int]]]]:
    """ Takes a string function name, returns the parser function """

    fn_map = {
        'sample': parse_sample_float,
        'monotonic': parse_monotonic_float,
        'seq': parse_seq_float,
    }

    return fn_map.get(fn, None)


def int_list_fn_finder(
    fn: str
) -> Optional[Callable[[Sequence[Any]], AOList[int]]]:
    """ Takes a string function name, returns the parser function """

    fn_map = {
        'sample': parse_sample_int,
        'monotonic': parse_monotonic_int,
        'seq': parse_seq_int,
    }

    return fn_map.get(fn, None)


def str_list_fn_finder(
    fn: str
) -> Optional[Callable[[Sequence[Any]], AOList[str]]]:
    """ Takes a string function name, returns the parser function """

    if fn == "sample":
        return parse_sample_str
    else:
        return None


# Distribution parsers


def parse_string_as_int(elem: Union[int, str]) -> Distribution[int]:
    try:
        parsed_elem: Distribution[int] = dist.IntConst(int(elem))
    except ValueError:
        raise DistParseError(f"Could not parse value '{elem}' as an integer.")

    return parsed_elem


def parse_list_in_dist(
    elems: Sequence[Any],
    fn_finder: Callable[[str], Optional[Callable[[Sequence[Any]], T]]],
) -> T:
    """ Parse a list while parsing a distribution.

    This could be a function or an actual list.
    """

    if len(elems) == 0:
        raise DistParseError("Encountered an empty list.")

    elif not isinstance(elems[0], str):
        raise DistParseError(
            f"Received a list but cannot take a list here: {elems}"
        )

    func = fn_finder(elems[0])

    if func is None:
        raise DistParseError(
            f"Recieved a list or invalid function: {elems}."
        )
    else:
        parsed_elem = func(elems[1:])

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
        parsed_elem = parse_list_in_dist(elem, int_dist_fn_finder)

    elif isinstance(elem, Distribution):
        parsed_elem = elem

    elif isinstance(elem, float):
        raise DistParseError(
            "This distribution/parameter doesn't support float values. "
            "Please convert '{elem}' to an integer."
        )

    else:
        raise DistParseError(
            "Recieved invalid input: {repr(elem)}."
        )

    if isinstance(parsed_elem, (dist.PMF, dist.Trunc)):
        return parsed_elem
    else:
        raise DistParseError(
            f"Could not parse '{elem}' as a discrete distribution."
        )


def parse_string_as_float(elem: Union[int, float, str]) -> DistributionIF:
    try:
        parsed_elem: DistributionIF = dist.FloatConst(float(elem))
    except ValueError:
        raise DistParseError(f"Could not parse value '{elem}' as a float.")
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
        parsed_elem = parse_list_in_dist(elem, float_dist_fn_finder)

    elif isinstance(elem, Distribution):
        parsed_elem = elem

    else:
        raise DistParseError(
            "Recieved invalid input: {repr(elem)}."
        )

    if isinstance(parsed_elem, (dist.PDF, dist.PMF, dist.Trunc)):
        return parsed_elem
    else:
        raise DistParseError(
            f"Could not parse '{elem}' as a continuous or "
            "discrete distribution."
        )


def parse_str_dist(elem: Any) -> Distribution[str]:
    """ Parses lists or literals to a categorical distribution. """

    if isinstance(elem, (int, float, str)):
        parsed_elem: Distribution[str] = dist.StrConst(str(elem))

    elif isinstance(elem, list):
        parsed_elem = parse_list_in_dist(elem, str_dist_fn_finder)

    elif isinstance(elem, Distribution):
        parsed_elem = elem

    else:
        raise DistParseError(
            "Recieved invalid input: {repr(elem)}."
        )

    if isinstance(parsed_elem, dist.Categorical):
        return parsed_elem

    else:
        raise DistParseError(
            f"Could not parse value '{elem}' as a categorical distribution."
        )


# List parsers
#
# These are types that return a list of distributions, rather than converting
# a list into a single distribution.

def parse_int_list(elems: Sequence[Any]) -> AOList[int]:
    """ Parse a list of integer distributions.

    Examples:
    >>> parse_int_list([1, 2, 3])
    AOList([IntConst(1), IntConst(2), IntConst(3)])

    >>> parse_int_list(["sample", 5, ["range", 1, 10]])
    Sample(5, Range(IntConst(1), IntConst(10)))

    >>> parse_int_list(["monotonic", 5, ["range", 1, 10]])
    Monotonic(5, Range(IntConst(1), IntConst(10)))

    >>> parse_int_list(["seq", 5, ["range", 1, 10]])
    Seq(5, Range(IntConst(1), IntConst(10)))
    """

    if len(elems) > 0 and isinstance(elems[0], str):
        func = int_list_fn_finder(elems[0])

        if func is not None:
            return func(elems[1:])

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

    >>> parse_float_list(["sample", 5, ["chi", 1, 1, 1]])
    Sample(5, Chi(FloatConst(1.0), FloatConst(1.0), FloatConst(1.0)))

    >>> parse_float_list(["monotonic", 5, ["chi", 1, 1, 1]])
    Monotonic(5, Chi(FloatConst(1.0), FloatConst(1.0), FloatConst(1.0)))

    >>> parse_float_list(["seq", 5, ["chi", 1, 1, 1]])
    Seq(5, Chi(FloatConst(1.0), FloatConst(1.0), FloatConst(1.0)))
    """

    if len(elems) > 0 and isinstance(elems[0], str):
        func = float_list_fn_finder(elems[0])

        if func is not None:
            return func(elems[1:])

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
    """ Parse a list of string distributions.

    >>> parse_str_list(["sample", 5, ["choose", ["one", "two", "three"]]])
    Sample(5, ChooseS([StrConst('one'), StrConst('two'), StrConst('three')]))

    >>> parse_str_list(["sample", 5, ["one", "two", "three"]])
    Traceback (most recent call last):
        ...
    DistParseError: Received a list or invalid function ...
    """

    if len(elems) > 0 and isinstance(elems[0], str):
        func = str_list_fn_finder(elems[0])

        if func is not None:
            return func(elems[1:])

    try:
        return AOList([parse_str_dist(e) for e in elems])
    except DistParseError as e:
        raise DistParseError(
            f"Error while parsing str list {repr(elems)}. {e.msg}"
        )


# Function argument parsers


def check_float_param(
    index: int,
    elems: Sequence[Any],
    name: str,
    func: str,
) -> DistributionIF:
    """ Does error handling and conversion of float parameters. """

    try:
        param = parse_float_dist(elems[index])
    except IndexError:
        raise DistParseError(
            f"Missing parameter <{name}> for function '{func}'. "
            f"You provided {elems}."
        )
    except ValueError:
        raise DistParseError(
            f"Could not parse parameter <{name}> value '{elems[index]}' "
            f"as a continuous or discrete distribution for function '{func}'."
        )
    return param


def check_int_param(
    index: int,
    elems: Sequence[Any],
    name: str,
    func: str
) -> Distribution[int]:
    """ Does error handling and conversion of int parameters. """

    try:
        param = parse_int_dist(elems[index])
    except IndexError:
        raise DistParseError(
            f"Missing parameter <{name}> for function '{func}'. "
            f"You provided {elems}."
        )
    except ValueError:
        raise DistParseError(
            f"Could not parse parameter <{name}> value '{elems[index]}' "
            f"as a discrete distribution for function '{func}'."
        )
    return param


def check_str_param(
    index: int,
    elems: Sequence[Any],
    name: str,
    func: str
) -> Distribution[str]:
    """ Does error handling and conversion of str parameters. """

    try:
        param = parse_str_dist(elems[index])
    except IndexError:
        raise DistParseError(
            f"Missing parameter <{name}> for function '{func}'. "
            f"You provided {elems}."
        )
    except ValueError:
        raise DistParseError(
            f"Could not parse parameter <{name}> value '{elems[index]}' "
            f"as a categorical distribution for function '{func}'."
        )
    return param


def check_true_int_param(
    index: int,
    elems: Sequence[Any],
    param: str,
    func: str,
) -> int:
    """ Parse an integer parameter that isn't allowed to be a distribution. """

    try:
        value = elems[index]
    except IndexError:
        raise DistParseError(
            f"Missing parameter <{param}> for function '{func}'. "
            f"You provided {elems}."
        )

    if isinstance(value, float):
        raise DistParseError(
            f"Parameter <{param}> value '{value}' in function '{func}' "
            "cannot be a float."
        )

    try:
        return int(value)
    except (ValueError, TypeError):
        raise DistParseError(
            f"Could not parse parameter <{param}> value '{elems[index]}' "
            f"as an integer in function '{func}'."
        )


def check_true_float_param(
    index: int,
    elems: Sequence[Any],
    param: str,
    func: str
) -> float:
    """ Parse a float parameter that isn't allowed to be a distribution. """

    try:
        value = elems[index]
    except IndexError:
        raise DistParseError(
            f"Missing parameter <{param}> for function '{func}'. "
            f"You provided {elems}."
        )

    try:
        return float(value)
    except (ValueError, TypeError):
        raise DistParseError(
            f"Could not parse parameter <{param}> value '{elems[index]}' "
            f"as a float in function '{func}'."
        )


def check_n_params(
    valid: Sequence[str],
    elems: Sequence[Any],
    func: str
) -> None:
    """ Checks if there are more parameters that was expected. """

    if len(elems) != len(valid):
        raise DistParseError(
            f"Function '{func}' expects {len(valid)} parameters: {valid}. "
            f"Received {len(elems)} parameters: {elems}."
        )
    return None


# Continuous distribution parsers


def parse_uniform(elems: Sequence[Any]) -> dist.Uniform:
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

    >>> parse_uniform(["1", [3, 2]])
    Traceback (most recent call last):
        ...
    DistParseError: Received a list but cannot take a list here: [3, 2]
    """

    str_func = "uniform"

    check_n_params(["min", "max"], elems, str_func)

    min_ = check_float_param(0, elems, "min", str_func)
    max_ = check_float_param(1, elems, "max", str_func)

    return dist.Uniform(min_, max_)


def parse_expuniform(elems: Sequence[Any]) -> dist.ExpUniform:
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

    str_func = "expuniform"

    check_n_params(["base", "min", "max"], elems, str_func)

    base = check_float_param(0, elems, "base", str_func)
    min_ = check_float_param(1, elems, "min", str_func)
    max_ = check_float_param(2, elems, "max", str_func)

    return dist.ExpUniform(base, min_, max_)


def parse_beta(elems: Sequence[Any]) -> dist.Beta:
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

    str_func = "beta"

    check_n_params(["a", "b"], elems, str_func)

    a = check_float_param(0, elems, "a", str_func)
    b = check_float_param(1, elems, "b", str_func)

    return dist.Beta(a, b)


def parse_chi(elems: Sequence[Any]) -> dist.Chi:
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

    str_func = "chi"

    check_n_params(["df", "loc", "scale"], elems, str_func)

    df = check_float_param(0, elems, "df", str_func)
    loc = check_float_param(1, elems, "loc", str_func)
    scale = check_float_param(2, elems, "scale", str_func)

    return dist.Chi(df, loc, scale)


def parse_gamma(elems: Sequence[Any]) -> dist.Gamma:
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

    str_func = "gamma"

    check_n_params(["a", "loc", "scale"], elems, str_func)

    a = check_float_param(0, elems, "a", str_func)
    loc = check_float_param(1, elems, "loc", str_func)
    scale = check_float_param(2, elems, "scale", str_func)

    return dist.Gamma(a, loc, scale)


def parse_normal(elems: Sequence[Any]) -> dist.Normal:
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

    str_func = "normal"

    check_n_params(["loc", "scale"], elems, str_func)

    loc = check_float_param(0, elems, "loc", str_func)
    scale = check_float_param(1, elems, "scale", str_func)

    return dist.Normal(loc, scale)


# Discrete distributions


def parse_range(elems: Sequence[Any]) -> dist.Range:
    """ Parse a list definition as a range distribution definition.

    Examples:

    >>> r = parse_range([1, 3])
    >>> r.min
    IntConst(1)

    >>> r.max
    IntConst(3)

    >>> r = parse_range([1, 2, 3])
    Traceback (most recent call last):
        ...
    DistParseError: ...
    """

    str_func = "range"

    check_n_params(["min", "max"], elems, str_func)

    min_ = check_int_param(0, elems, "min", str_func)
    max_ = check_int_param(1, elems, "max", str_func)

    # if min_ >= max_:
    #     raise DistParseError(f"The <min> bound must be smaller than <max>")

    return dist.Range(min_, max_)


def parse_binomial(elems: Sequence[Any]) -> dist.Binomial:
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

    str_func = "binomial"

    check_n_params(["n", "p", "loc"], elems, str_func)

    n = check_int_param(0, elems, "n", str_func)
    p = check_float_param(1, elems, "p", str_func)

    if not isinstance(p, dist.PDF):
        raise DistParseError(
            "Parameter <p> to binomial must be a float between 0 and 1. "
            "You provided '{elems[1]}'."
        )

    loc = check_int_param(2, elems, "loc", str_func)

    return dist.Binomial(n, p, loc)


def parse_neg_binomial(elems: Sequence[Any]) -> dist.NegBinomial:
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

    str_func = "nbinomial"

    check_n_params(["n", "p", "loc"], elems, str_func)

    n = check_int_param(0, elems, "n", str_func)
    p = check_float_param(1, elems, "p", str_func)

    if not isinstance(p, dist.PDF):
        raise DistParseError(
            "Parameter <p> to nbinomial must be a float between 0 and 1. "
            "You provided '{elems[1]}'."
        )

    loc = check_int_param(2, elems, "loc", str_func)

    return dist.NegBinomial(n, p, loc)


def parse_hypergeom(elems: Sequence[Any]) -> dist.HyperGeometric:
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
    str_func = "hypergeom"

    check_n_params(
        ["n_success", "n", "population", "loc"],
        elems,
        str_func
    )

    n_success = check_int_param(0, elems, "n_success", str_func)
    n = check_int_param(1, elems, "n", str_func)
    population = check_int_param(2, elems, "population", str_func)
    loc = check_int_param(3, elems, "loc", str_func)

    return dist.HyperGeometric(n_success, n, population, loc)


# Distribution transformers.


def parse_trunc_int(elems: Sequence[Any]) -> dist.Trunc[int]:
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

    str_func = "trunc"

    check_n_params(["min", "max", "dist"], elems, str_func)

    min_ = check_true_int_param(0, elems, "min", str_func)
    max_ = check_true_int_param(1, elems, "max", str_func)
    this_dist = check_int_param(2, elems, "dist", str_func)

    return dist.Trunc(min_, max_, this_dist)


def parse_trunc_float(elems: Sequence[Any]) -> dist.Trunc[Union[float, int]]:
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

    str_func = "trunc"

    check_n_params(["min", "max", "dist"], elems, str_func)

    min_ = check_true_float_param(0, elems, "min", str_func)
    max_ = check_true_float_param(1, elems, "max", str_func)

    # I'm being a bit cheeky with type hints here.
    # Return from check float param is DistributionIF
    # The type system isn't expressive enough to destructure to
    # Trunc[Union[int, float]].
    this_dist: Distribution[Any] = check_float_param(2, elems, "dist", str_func)
    return dist.Trunc(min_, max_, this_dist)


def parse_choose_int(elems: Sequence[Any]) -> dist.ChooseI:
    """ Parse choose from a list representation.

    Examples:
    >>> c = parse_choose_int([[1, 2, ['range', 5, 10]]])
    >>> c.choices
    [IntConst(1), IntConst(2), Range(IntConst(5), IntConst(10))]
    """

    check_n_params(["dist"], elems, "choose")

    parsed_elems = parse_int_list(elems[0])

    return dist.ChooseI(parsed_elems.elems)


def parse_choose_float(elems: Sequence[Any]) -> dist.ChooseF:
    """ Parse choose from a list representation.

    Examples:
    >>> c = parse_choose_float([[2, ['uniform', 5, 10]]])
    >>> c.choices
    [FloatConst(2.0), Uniform(FloatConst(5.0), FloatConst(10.0))]
    """

    check_n_params(["dist"], elems, "choose")

    parsed_elems = parse_float_list(elems[0])

    return dist.ChooseF(parsed_elems.elems)


def parse_choose_str(elems: Sequence[Any]) -> dist.ChooseS:
    """ Parse choose from a list representation.

    Examples:
    >>> c = parse_choose_str([[2, "two"]])
    >>> c.choices
    [StrConst('2'), StrConst('two')]
    """

    check_n_params(["dist"], elems, "choose")

    parsed_elems = parse_str_list(elems[0])

    return dist.ChooseS(parsed_elems.elems)


# List function parsers


def parse_sample_float(elems: Sequence[Any]) -> AOList[Union[int, float]]:
    """ Choose n floats from a distribution.

    Examples:

    >>> c = parse_sample_float([4, ["uniform", 0, 1]])
    >>> c
    Sample(4, Uniform(FloatConst(0.0), FloatConst(1.0)))

    >>> c = parse_sample_float([4, ["range", 0, 10]])
    >>> c
    Sample(4, Range(IntConst(0), IntConst(10)))
    """

    str_func = "sample"

    check_n_params(["n", "dist"], elems, str_func)

    n = check_true_int_param(0, elems, "n", str_func)
    this_dist: Distribution = check_float_param(1, elems, "dist", str_func)
    return dist.Sample(n, this_dist)


def parse_sample_int(elems: Sequence[Any]) -> AOList[int]:
    """ Choose n floats from a distribution.

    Examples:

    >>> c = parse_sample_int([4, ["range", 0, 10]])
    >>> c
    Sample(4, Range(IntConst(0), IntConst(10)))

    >>> c = parse_sample_int([4, ["uniform", 0, 1]])
    Traceback (most recent call last):
        ...
    DistParseError: Recieved a list or invalid function ['uniform', 0, 1].
    """

    str_func = "sample"

    check_n_params(["n", "dist"], elems, str_func)

    n = check_true_int_param(0, elems, "n", str_func)
    this_dist = check_int_param(1, elems, "dist", str_func)
    return dist.Sample(n, this_dist)


def parse_sample_str(elems: Sequence[Any]) -> AOList[str]:
    """ Choose n floats from a distribution.

    Examples:

    >>> c = parse_sample_str([4, ["choose", ["one", "two"]]])
    >>> c
    Sample(4, ChooseS([StrConst('one'), StrConst('two')]))
    """

    str_func = "sample"

    check_n_params(["n", "dist"], elems, str_func)

    n = check_true_int_param(0, elems, "n", str_func)
    this_dist = check_str_param(1, elems, "dist", str_func)

    return dist.Sample(n, this_dist)


def parse_monotonic_int(elems: Sequence[Any]) -> AOList[int]:
    """ Choose n ints from distribution and sort them in increasing order.

    Examples:

    >>> parse_monotonic_int([4, ["binomial", 4, 0.2, 0]])
    Monotonic(4, Binomial(IntConst(4), FloatConst(0.2), IntConst(0)))
    """

    str_func = "monotonic"

    check_n_params(["n", "dist"], elems, str_func)

    n = check_true_int_param(0, elems, "n", str_func)
    this_dist = check_int_param(1, elems, "dist", str_func)
    return dist.Monotonic(n, this_dist)


def parse_monotonic_float(elems: Sequence[Any]) -> AOList[Union[int, float]]:
    """ Choose n floats from distribution and sort them in increasing order.

    Examples:

    >>> parse_monotonic_float([4, ["uniform", 0, 1]])
    Monotonic(4, Uniform(FloatConst(0.0), FloatConst(1.0)))
    """

    str_func = "monotonic"

    check_n_params(["n", "dist"], elems, str_func)

    n = check_true_int_param(0, elems, "n", str_func)
    this_dist: Distribution = check_float_param(1, elems, "dist", str_func)
    return dist.Monotonic(n, this_dist)


def parse_seq_int(elems: Sequence[Any]) -> AOList[int]:
    """ Choose n quantiles from distribution.

    Examples:

    >>> parse_seq_int([4, ["binomial", 10, 0.2, 0]])
    Seq(4, Binomial(IntConst(10), FloatConst(0.2), IntConst(0)))
    """

    str_func = "seq"

    check_n_params(["n", "dist"], elems, str_func)

    n = check_true_int_param(0, elems, "n", str_func)
    this_dist = check_int_param(1, elems, "dist", str_func)
    return dist.Seq(n, this_dist)


def parse_seq_float(elems: Sequence[Any]) -> AOList[Union[int, float]]:
    """ Choose n quantiles from distribution.

    Examples:

    >>> parse_seq_float([4, ["uniform", 0, 1]])
    Seq(4, Uniform(FloatConst(0.0), FloatConst(1.0)))
    """

    str_func = "seq"

    check_n_params(["n", "dist"], elems, str_func)

    n = check_true_int_param(0, elems, "n", str_func)
    this_dist: Distribution = check_float_param(1, elems, "dist", str_func)
    return dist.Seq(n, this_dist)
