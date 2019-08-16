"""
"""

import re
from typing import cast
from typing import List
from typing import Sequence
from typing import Union


from augustus_optimiser.errors import DistParseError

WHITESPACE = re.compile(r"\s+")


class Distribution(object):
    """ Abstract base class for probability distribution functions.

    Not intedended to be used directly.
    Subclasses should implement these methods.
    """

    @classmethod
    def parse(cls, string: str) -> "Distribution":
        """ Parse a string as a function definition. """

        split_string = WHITESPACE.split(string.strip())
        return cls.parse_list(cast(List[Union[str, int, float]], split_string))

    @classmethod
    def parse_list(
        cls,
        li: Sequence[Union[str, int, float]]
    ) -> "Distribution":
        raise NotImplementedError("Subclasses must implement this method.")

    @staticmethod
    def usage() -> str:
        raise NotImplementedError("Subclasses must implement this method.")


class PDF(Distribution):
    """ Abstract base class from continuous distributions. """

    def random(self) -> float:
        raise NotImplementedError("Subclasses must implement this method.")

    def many_random(self, n: int) -> List[float]:
        raise NotImplementedError("Subclasses must implement this method.")


class PMF(Distribution):
    """ Abstract base class for discrete distributions. """

    def random(self) -> int:
        raise NotImplementedError("Subclasses must implement this method.")

    def many_random(self, n: int) -> List[int]:
        raise NotImplementedError("Subclasses must implement this method.")


class Categorical(Distribution):
    """ Abstract base class for selecting one of many choices. """

    def random(self) -> str:
        raise NotImplementedError("Subclasses must implement this method.")

    def many_random(self, n: int) -> List[str]:
        raise NotImplementedError("Subclasses must implement this method.")


class Uniform(PDF):

    def __init__(self, min: float, max: float) -> None:
        self.min = min
        self.max = max
        return

    @staticmethod
    def usage() -> str:
        return "uniform <min> <max>"

    @classmethod
    def parse_list(cls, li: Sequence[Union[str, int, float]]) -> "Uniform":
        """ Parse a list definition as a uniform distribution definition.

        Examples:

        >>> u = Uniform.parse_list(["uniform", 1, 3.5])
        >>> u.min
        1.0
        >>> u.max
        3.5
        >>> Uniform.parse_list(["", 1, 2])
        Traceback (most recent call last):
            ...
        DistParseError: Invalid function name for uniform: ''.
        """

        fn_name = li[0]
        if fn_name != "uniform":
            raise DistParseError(
                f"Invalid function name for uniform: '{fn_name}'",
            )

        try:
            min_ = float(li[1])
        except IndexError:
            raise DistParseError(
                "Missing min value for pdf.",
            )
        except ValueError:
            raise DistParseError(
                f"Could not parse min value as float: {li[1]}.",
            )

        try:
            max_ = float(li[2])
        except IndexError:
            raise DistParseError(
                "Missing max value for pdf.",
            )
        except ValueError:
            raise DistParseError(
                f"Could not parse max value as float: {li[2]}.",
            )

        if min_ >= max_:
            raise DistParseError(f"The <min> bound must be smaller than <max>")

        return cls(min_, max_)

    def random(self) -> float:
        """ Return a random value between min and max.

        Examples:
        >>> for _ in range(100):
        ...     x = Uniform(0, 5).random()
        ...     assert isinstance(x, float)
        ...     assert 0 <= x <= 5
        """

        import random
        return random.uniform(self.min, self.max)

    def many_random(self, n: int) -> List[float]:
        """ Return n random values between min and max.

        Examples:
        >>> x = Uniform(0, 5).many_random(100)
        >>> assert isinstance(x, list)
        >>> assert all(map(lambda y: isinstance(y, float), x))
        >>> assert all(map(lambda y: 0 <= y <= 5, x))
        """

        import random
        return [random.uniform(self.min, self.max) for _ in range(n)]


class Range(PMF):

    def __init__(self, min: int, max: int) -> None:
        self.min = min
        self.max = max
        return

    @staticmethod
    def usage() -> str:
        return "range <min> <max>"

    @classmethod
    def parse_list(cls, li: Sequence[Union[str, int, float]]) -> "Range":
        """ Parse a list definition as a range distribution definition.

        Examples:

        >>> r = Range.parse_list(["range", 1, 3])
        >>> r.min
        1
        >>> r.max
        3
        >>> Range.parse_list(["hi", 1, 3])
        Traceback (most recent call last):
            ...
        DistParseError: Invalid function name for uniform: 'hi'.
        >>> Range.parse_list(["range", 1, 3.5])
        Traceback (most recent call last):
            ...
        DistParseError: Max bound is a float: 3.5. Expected an integer.
        """

        fn_name = li[0]
        if fn_name != "range":
            raise DistParseError(
                f"Invalid function name for range: '{fn_name}'",
            )

        try:
            min_ = int(li[1])
            assert not isinstance(li[1], float)
        except IndexError:
            raise DistParseError(
                "Missing min value for range.",
            )
        except AssertionError:
            raise DistParseError(f"Min bound is a float: {li[1]}. "
                                 "Expected an integer.")
        except ValueError:
            raise DistParseError(
                f"Could not parse min value as integer: {li[1]}.",
            )

        try:
            max_ = int(li[2])
            assert not isinstance(li[2], float)
        except IndexError:
            raise DistParseError(
                "Missing max value for range.",
            )
        except AssertionError:
            raise DistParseError(f"Max bound is a float: {li[2]}. "
                                 "Expected an integer.")
        except ValueError:
            raise DistParseError(
                f"Could not parse max value as integer: {li[2]}.",
            )

        if min_ >= max_:
            raise DistParseError(f"The <min> bound must be smaller than <max>")

        return cls(min_, max_)

    def random(self) -> int:
        """ Return a random integer between min and max.

        Examples:
        >>> for _ in range(100):
        ...     x = Range(0, 5).random()
        ...     assert isinstance(x, int)
        ...     assert 0 <= x <= 5
        """

        import random
        return random.randint(self.min, self.max)

    def many_random(self, n: int) -> List[int]:
        """ Return n random values between min and max.

        Examples:
        >>> x = Range(0, 5).many_random(100)
        >>> assert all(map(lambda y: isinstance(y, int), x))
        >>> assert all(map(lambda y: 0 <= y <= 5, x))
        """

        import random
        return [random.randint(self.min, self.max) for _ in range(n)]


class Choose(Categorical):

    def __init__(self, choices: Sequence[str]) -> None:
        self.choices = list(choices)
        return

    @staticmethod
    def usage() -> str:
        return "choose option1 option2 ..."

    @classmethod
    def parse_list(cls, li: Sequence[Union[str, int, float]]) -> "Choose":
        """ Parse choose from a list representation.

        Examples:
        >>> c = Choose.parse_list(["choose", "one", "two", 3])
        >>> c.choices
        ['one', 'two', '3']
        >>> Choose.parse_list(["", 1, 2])
        Traceback (most recent call last):
            ...
        DistParseError: Invalid function name for choose: ''.
        >>> Choose.parse_list(["choose"])
        Traceback (most recent call last):
            ...
        DistParseError: No options specified for choose to use..
        """

        fn_name = li[0]
        if fn_name != "choose":
            raise DistParseError(
                f"Invalid function name for choose: '{fn_name}'",
            )

        choices = li[1:]
        if len(choices) == 0:
            raise DistParseError("No options specified for choose to use.")

        return cls(list(map(str, choices)))

    def random(self) -> str:
        """ Returns a random choice.

        >>> x = Choose(["one", "two", "three"]).random()
        >>> assert x in ["one", "two", "three"]
        """

        import random
        return random.choice(self.choices)

    def many_random(self, n: int) -> List[str]:
        """ Returns n random choices.

        Examples:
        >>> x = Choose(["one", "two", "three"]).many_random(5)
        >>> assert all(map(lambda y: y in ["one", "two", "three"], x))
        """

        import random
        return random.choices(self.choices, k=n)
