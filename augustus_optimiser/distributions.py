"""
"""

from typing import List
from typing import Sequence
from typing import Union
from typing import Generic
from typing import TypeVar

from augustus_optimiser.errors import DistValueError


# Genetic type var for baseclassing Distribution.
T = TypeVar("T", str, int, float)


class Distribution(Generic[T]):
    """ Abstract base class for probability distribution functions.

    Not intedended to be used directly.
    Subclasses should implement these methods.
    """

    def get(self) -> T:
        raise NotImplementedError("Subclasses must implement this method.")

    def get_many(self, n: int) -> List[T]:
        return [self.get() for _ in range(n)]


class FloatConst(Distribution[float]):
    """ A constant value. """

    def __init__(self, value: float) -> None:
        self.value = value
        return

    def get(self) -> float:
        return self.value


class IntConst(Distribution[int]):

    def __init__(self, value: int) -> None:
        self.value = value
        return

    def get(self) -> int:
        return self.value


class StrConst(Distribution[str]):

    def __init__(self, value: str) -> None:
        self.value = value
        return

    def get(self) -> str:
        return self.value


PDFUnion = Union[float, Distribution[float]]


class Uniform(Distribution[float]):

    def __init__(self, min: PDFUnion, max: PDFUnion) -> None:
        if isinstance(min, float):
            self.min: Distribution[float] = FloatConst(min)
        else:
            self.min = min

        if isinstance(max, float):
            self.max: Distribution[float] = FloatConst(max)
        else:
            self.max = max
        return

    def get(self) -> float:
        """ Return a random value between min and max.

        Examples:
        >>> for _ in range(100):
        ...     x = Uniform(0, 5).get()
        ...     assert isinstance(x, float)
        ...     assert 0 <= x <= 5
        """

        import random
        return random.uniform(self.min.get(), self.max.get())


PMFUnion = Union[int, Distribution[int]]


class Range(Distribution[int]):

    def __init__(self, min: PMFUnion, max: PMFUnion) -> None:
        if isinstance(min, int):
            self.min: Distribution[int] = IntConst(min)
        else:
            self.min = min

        if isinstance(max, int):
            self.max: Distribution[int] = IntConst(max)
        else:
            self.max = max
        return

    def get(self) -> int:
        """ Return a random integer between min and max.

        Examples:
        >>> for _ in range(100):
        ...     x = Range(0, 5).random()
        ...     assert isinstance(x, int)
        ...     assert 0 <= x <= 5
        """

        import random
        return random.randint(self.min.get(), self.max.get())


class Choose(Distribution[T]):

    def __init__(self, choices: Sequence[T]) -> None:
        self.choices: List[T] = list(choices)

        first_type = type(self.choices[0])
        if not all(isinstance(t, first_type) for t in self.choices):
            raise DistValueError(
                "All members of categorical distributions must be the same "
                "type."
            )
        return

    def get(self) -> T:
        """ Returns a random choice.

        >>> x = Choose(["one", "two", "three"]).random()
        >>> assert x in ["one", "two", "three"]
        """

        import random
        return random.choice(self.choices)

    def get_many(self, n: int) -> List[T]:
        """ Returns n random choices.

        Examples:
        >>> x = Choose(["one", "two", "three"]).many_random(5)
        >>> assert all(map(lambda y: y in ["one", "two", "three"], x))
        """

        import random
        return random.choices(self.choices, k=n)


# monotonic 3 range 1 5
# 2 4 5

# sample 3 range 1 5
# 3 1 5

# seq 5 uniform 0.5 0.9
# 0.5 0.6 0.7 0.8 0.9
