from typing import Any
from typing import Sequence
from typing import TypeVar


from augustus_optimiser.distributions import Distribution
from augustus_optimiser.distributions import FloatConst
from augustus_optimiser.distributions import IntConst
from augustus_optimiser.distributions import StrConst
from augustus_optimiser.distributions import Uniform
from augustus_optimiser.distributions import Range
from augustus_optimiser.distributions import Choose
from augustus_optimiser.errors import DistParseError


T = TypeVar("T", str, int, float)

# STR_TO_PARSER = {
#     "uniform": parse_uniform,
# }


def parse(elems: Sequence[Any]) -> Distribution:
    return elems[0]


def parse_float_const(elems: Sequence[Any]) -> FloatConst:
    if len(elems) == 0:
        raise DistParseError("Missing value, expected float.")
    elif len(elems) > 1:
        raise DistParseError("Received multiple values but "
                             "expected a single float.")

    if isinstance(elems[0], list):
        elem = parse(elems[0])
    else:
        elem = elems[0]

    if isinstance(elem, FloatConst):
        return FloatConst(elem.value)

    if isinstance(elem, (int, float, str)):
        try:
            elem = float(elem)
        except ValueError:
            raise DistParseError(f"Could not parse value 'elems[0]' as float.")

    else:
        raise DistParseError(f"Could not parse value 'elems[0]' as float.")

    return FloatConst(elem)


def parse_int_const(elems: Sequence[Any]) -> IntConst:
    if len(elems) == 0:
        raise DistParseError("Missing value, expected int.")
    elif len(elems) > 1:
        raise DistParseError("Received multiple values but "
                             "expected a single int.")

    if isinstance(elems[0], list):
        elem = parse(elems[0])
    else:
        elem = elems[0]

    if isinstance(elem, IntConst):
        return IntConst(elem.value)

    if isinstance(elem, (int, str)):
        try:
            elem = int(elem)
        except (ValueError, TypeError):
            raise DistParseError(f"Could not parse value 'elems[0]' as int.")

    else:
        raise DistParseError(f"Could not parse value 'elems[0]' as int.")

    return IntConst(elem)


def parse_str_const(elems: Sequence[Any]) -> StrConst:
    if len(elems) == 0:
        raise DistParseError("Missing value, expected str.")
    elif len(elems) > 1:
        raise DistParseError("Received multiple values but "
                             "expected a str float.")

    if isinstance(elems[0], list):
        elem = parse(elems[0])
    else:
        elem = elems[0]

    try:
        fl = str(elem)
    except (ValueError, TypeError):
        raise DistParseError(f"Could not parse value 'elems[0]' as str.")

    return StrConst(fl)


def parse_uniform(elems: Sequence[Any]) -> Uniform:
    """ Parse a list definition as a uniform distribution definition.

    Examples:

    >>> u, _ = Uniform.parse_list([1, 3.5])
    >>> u.min
    1.0
    >>> u.max
    3.5
    >>> Uniform.parse_list(["one", 2])
    Traceback (most recent call last):
        ...
    DistParseError: Could not parse min value 'one' as float.
    """

    try:
        min_ = float(elems[0])
    except IndexError:
        raise DistParseError(
            "Missing <min> parameter for pdf.",
        )
    except ValueError:
        raise DistParseError(
            f"Could not parse <min> parameter '{elems[0]}' as float.",
        )

    try:
        max_ = float(elems[1])
    except IndexError:
        raise DistParseError(
            "Missing <max> parameter for pdf.",
        )
    except ValueError:
        raise DistParseError(
            f"Could not parse <max> parameter '{elems[1]}' as float.",
        )

    if min_ >= max_:
        raise DistParseError(f"The <min> bound must be smaller than <max>")

    return Uniform(min_, max_)


def parse_range(elems: Sequence[Any]) -> Range:
    """ Parse a list definition as a range distribution definition.

    Examples:

    >>> r, _ = Range.parse_list([1, 3])
    >>> r.min
    1
    >>> r.max
    3
    >>> r, l = Range.parse_list([1, 2, 3])
    >>> r.max
    2
    >>> l
    [3]
    >>> Range.parse_list(["range", 1])
    Traceback (most recent call last):
        ...
    DistParseError: Could not parse min parameter 'range' as integer.
    """

    try:
        min_ = int(elems[0])
        assert not isinstance(elems[0], float)
    except IndexError:
        raise DistParseError(
            "Missing min parameter for range.",
        )
    except (ValueError, AssertionError):
        raise DistParseError(
            f"Could not parse min parameter '{elems[0]}' as integer.",
        )

    try:
        max_ = int(elems[1])
        assert not isinstance(elems[1], float)
    except IndexError:
        raise DistParseError(
            "Missing max parameter for range.",
        )
    except (AssertionError, ValueError):
        raise DistParseError(
            f"Could not parse max parameter '{elems[1]}' as integer.",
        )

    if min_ >= max_:
        raise DistParseError(f"The <min> bound must be smaller than <max>")

    return Range(min_, max_)


def parse_choose(elems: Sequence[Any]) -> Choose[T]:
    """ Parse choose from a list representation.

    Examples:
    >>> c, r = Choose.parse_list(["one", "two", 3])
    >>> c.choices
    ['one', 'two', '3']
    >>> r
    []
    """

    if len(elems) == 0:
        raise DistParseError("No options specified for choose to use.")

    assert False  # Figure out what type it should be.

    return Choose(elems)
