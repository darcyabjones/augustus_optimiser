"""
"""

from typing import List
from typing import Set
from typing import Union
from typing import Sequence
from typing import Generic
from typing import TypeVar
from typing import Iterable


import numpy as np
from scipy import stats


from ao.errors import DistDomainError, DistValueError


# Genetic type var for baseclassing Distribution.
T = TypeVar("T", str, int, float, Union[int, float])
U = TypeVar('U')
V = TypeVar("V", int, float, Union[int, float])


# Base types for distributions


class Distribution(Generic[T]):
    """ Abstract base class for probability distribution functions.

    Not intedended to be used directly.
    Subclasses should implement these methods.
    """

    def get(self) -> T:
        raise NotImplementedError("Subclasses must implement this method.")

    def get_many(self, n: int) -> List[T]:
        return [self.get() for _ in range(n)]

    def seq(self, n: int) -> List[T]:
        """ Should return a non-redundant grid across the distribution """

        raise NotImplementedError("Subclasses must implement this method.")


# We can use integer distributions in float distributions.
DistributionIF = Union[Distribution[int], Distribution[float]]


class PDF(Distribution[float]):
    """ Concrete type for a distribution of floats. """
    pass


class PMF(Distribution[int]):
    """ Concrete type for a distribution of ints. """
    pass


class Categorical(Distribution[T]):

    def seq(self, n: int) -> List[T]:
        raise DistValueError(
            "Categorical distributions do not support the seq operation."
        )


# Atomic values
# We use these to allow using .get() everywhere.


class FloatConst(PDF):
    """ A constant float value. """

    def __init__(self, value: Union[float]) -> None:
        self.value = float(value)
        return

    def get(self) -> float:
        return self.value

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({repr(self.value)})"

    def seq(self, n: int) -> List[float]:
        raise DistValueError("Atomic values do not support the seq operation.")


class IntConst(PMF):
    """ A constant integer value. """

    def __init__(self, value: int) -> None:
        self.value = value
        return

    def get(self) -> int:
        return self.value

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({repr(self.value)})"

    def seq(self, n: int) -> List[int]:
        raise DistValueError("Atomic values do not support the seq operation.")


class StrConst(Categorical[str]):
    """ A constant string value. """

    def __init__(self, value: str) -> None:
        self.value = value
        return

    def get(self) -> str:
        return self.value

    def __repr__(self):
        return f"{self.__class__.__name__}({repr(self.value)})"

    def seq(self, n: int) -> List[str]:
        raise DistValueError("Atomic values do not support the seq operation.")


# Distributions


class Uniform(PDF):

    def __init__(self, min: DistributionIF, max: DistributionIF) -> None:
        """ Take a random float between min and max """

        if isinstance(min, float):
            self.min: DistributionIF = FloatConst(min)
        elif isinstance(min, int):
            self.min = FloatConst(float(min))
        else:
            self.min = min

        if isinstance(max, float):
            self.max: DistributionIF = FloatConst(max)
        elif isinstance(max, int):
            self.max = FloatConst(float(max))
        else:
            self.max = max
        return

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({repr(self.min)}, {repr(self.max)})"

    def get(self) -> float:
        """ Return a random value between min and max.

        Examples:
        >>> for _ in range(100):
        ...     x = Uniform(0.0, 5.0).get()
        ...     assert isinstance(x, float)
        ...     assert 0 <= x <= 5
        """
        min_ = self.min.get()
        max_ = self.max.get()
        rv = stats.uniform.rvs(
            loc=min_,
            scale=(max_ - min_),
            size=None,
        )

        return float(rv)

    def get_many(self, n: int) -> List[float]:
        """ Return n random values between min and max.

        Examples:
        >>> for x in Uniform(0.0, 5.0).get_many(100):
        ...     assert isinstance(x, float)
        ...     assert 0 <= x <= 5
        """
        min_ = self.min.get()
        max_ = self.max.get()
        rv: List[float] = stats.uniform.rvs(
            loc=min_,
            scale=(max_ - min_),
            size=n,
        ).tolist()
        return rv

    def seq(self, n: int) -> List[float]:
        """ Return n quantiles between min and max.

        Examples:
        >>> for x in Uniform(0, 5).seq(100):
        ...     assert isinstance(x, float)
        ...     assert 0 <= x <= 5
        """

        min_ = self.min.get()
        max_ = self.max.get()
        return np.linspace(min_, max_, num=n).tolist()


class ExpUniform(PDF):

    def __init__(
        self,
        base: DistributionIF,
        min: DistributionIF,
        max: DistributionIF
    ) -> None:
        """ Take a random float between min and max on exponential scale. """

        if isinstance(base, float):
            self.base: DistributionIF = FloatConst(base)
        elif isinstance(base, int):
            self.base = FloatConst(float(base))
        else:
            self.base = base

        if isinstance(min, float):
            self.min: DistributionIF = FloatConst(min)
        elif isinstance(min, int):
            self.min = FloatConst(float(min))
        else:
            self.min = min

        if isinstance(max, float):
            self.max: DistributionIF = FloatConst(max)
        elif isinstance(max, int):
            self.max = FloatConst(float(max))
        else:
            self.max = max
        return

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({repr(self.min)}, {repr(self.max)})"

    def get(self) -> float:
        """ Return a random value on an exponential scale between min and max.

        Examples:
        >>> for _ in range(100):
        ...     x = ExpUniform(10, 0.0, 6.0).get()
        ...     assert isinstance(x, float)
        ...     assert 1 <= x <= 1e6
        """
        base = self.base.get()
        min_ = self.min.get()
        max_ = self.max.get()
        rv = stats.uniform.rvs(
            loc=min_,
            scale=(max_ - min_),
            size=None,
        )

        return float(base ** rv)

    def get_many(self, n: int) -> List[float]:
        """ Return n random values between min and max on exponential scale.

        Examples:
        >>> for x in ExpUniform(10, 0.0, 6.0).get_many(100):
        ...     assert isinstance(x, float)
        ...     assert 1 <= x <= 1e6
        """
        base = self.base.get()
        min_ = self.min.get()
        max_ = self.max.get()
        rv = stats.uniform.rvs(
            loc=min_,
            scale=(max_ - min_),
            size=n,
        )

        return (base ** rv).tolist()

    def seq(self, n: int) -> List[float]:
        """ Return n quantiles between min and max on exponential scale.

        Examples:
        >>> for x in ExpUniform(10, 0, 5).seq(100):
        ...     assert isinstance(x, float)
        ...     assert 1 <= x <= 1e5
        """

        base = self.base.get()
        min_ = self.min.get()
        max_ = self.max.get()
        return (base ** np.linspace(min_, max_, num=n)).tolist()


class Beta(PDF):

    def __init__(self, a: DistributionIF, b: DistributionIF) -> None:
        """ Take a random float from the beta distribution. """

        if isinstance(a, float):
            self.a: DistributionIF = FloatConst(a)
        elif isinstance(a, int):
            self.a = FloatConst(float(a))
        else:
            self.a = a

        if isinstance(b, float):
            self.b: DistributionIF = FloatConst(b)
        elif isinstance(b, int):
            self.b = FloatConst(float(b))
        else:
            self.b = b
        return

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({repr(self.a)}, {repr(self.b)})"

    def get(self) -> float:
        """ Return a random value from beta distribution.

        Examples:
        >>> for _ in range(100):
        ...     x = Beta(2.31, 0.627).get()
        ...     assert isinstance(x, float)
        ...     assert 0 <= x <= 1
        """
        a = self.a.get()
        b = self.b.get()
        rv = stats.beta.rvs(
            a=a,
            b=b,
            size=None,
        )

        return float(rv)

    def get_many(self, n: int) -> List[float]:
        """ Return n random values from the beta distribution.

        Examples:
        >>> for x in Beta(2.31, 0.627).get_many(100):
        ...     assert isinstance(x, float)
        ...     assert 0 <= x <= 1
        """
        a = self.a.get()
        b = self.b.get()
        rv = stats.beta.rvs(
            a=a,
            b=b,
            size=n,
        )

        return rv.tolist()

    def seq(self, n: int) -> List[float]:
        """ Return n quantiles from the beta distribution.

        Examples:
        >>> for x in Beta(5, 1).seq(100):
        ...     assert isinstance(x, float)
        ...     assert 0 <= x <= 1
        """

        a = self.a.get()
        b = self.b.get()

        quantiles = np.linspace(0.01, 0.99, num=n)
        return (stats.beta
                .ppf(quantiles, a, b)
                .astype(float)
                .tolist())


class Chi(PDF):

    def __init__(
        self,
        df: DistributionIF,
        loc: DistributionIF,
        scale: DistributionIF
    ) -> None:
        """ Take a random float from the chi distribution. """

        if isinstance(df, float):
            self.df: DistributionIF = FloatConst(df)
        elif isinstance(df, int):
            self.df = FloatConst(float(df))
        else:
            self.df = df

        if isinstance(loc, float):
            self.loc: DistributionIF = FloatConst(loc)
        elif isinstance(loc, int):
            self.loc = FloatConst(float(loc))
        else:
            self.loc = loc

        if isinstance(scale, float):
            self.scale: DistributionIF = FloatConst(scale)
        elif isinstance(loc, int):
            self.scale = FloatConst(float(scale))
        else:
            self.scale = scale
        return

    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}"
                f"({repr(self.df)}, {repr(self.loc)}, {repr(self.scale)})")

    def get(self) -> float:
        """ Return a random value from chi distribution.

        Examples:
        >>> for _ in range(100):
        ...     x = Chi(1, 0, 0.627).get()
        ...     assert isinstance(x, float)
        ...     assert 0 <= x
        """
        df = self.df.get()
        loc = self.loc.get()
        scale = self.scale.get()
        rv = stats.chi.rvs(
            df=df,
            loc=loc,
            scale=scale,
            size=None,
        )

        return float(rv)

    def get_many(self, n: int) -> List[float]:
        """ Return n random values from the chi distribution.

        Examples:
        >>> for x in Chi(1, 1, 100).get_many(100):
        ...     assert isinstance(x, float)
        ...     assert 1 <= x
        """
        df = self.df.get()
        loc = self.loc.get()
        scale = self.scale.get()
        rv = stats.chi.rvs(
            df=df,
            loc=loc,
            scale=scale,
            size=n,
        )

        return rv.tolist()

    def seq(self, n: int) -> List[float]:
        """ Return n quantiles from the chi distribution.

        Examples:
        >>> for x in Chi(1, 5, 1).seq(100):
        ...     assert isinstance(x, float)
        ...     assert 5 <= x
        """

        df = self.df.get()
        loc = self.loc.get()
        scale = self.scale.get()

        quantiles = np.linspace(0.01, 0.99, num=n)
        return (stats.chi
                .ppf(quantiles, df=df, loc=loc, scale=scale)
                .astype(float)
                .tolist())


class Gamma(PDF):

    def __init__(
        self,
        a: DistributionIF,
        loc: DistributionIF,
        scale: DistributionIF
    ) -> None:
        """ Take a random float from the gamma distribution. """

        if isinstance(a, float):
            self.a: DistributionIF = FloatConst(a)
        elif isinstance(a, int):
            self.a = FloatConst(float(a))
        else:
            self.a = a

        if isinstance(loc, float):
            self.loc: DistributionIF = FloatConst(loc)
        elif isinstance(loc, int):
            self.loc = FloatConst(float(loc))
        else:
            self.loc = loc

        if isinstance(scale, float):
            self.scale: DistributionIF = FloatConst(scale)
        elif isinstance(loc, int):
            self.scale = FloatConst(float(scale))
        else:
            self.scale = scale
        return

    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}"
                f"({repr(self.a)}, {repr(self.loc)}, {repr(self.scale)})")

    def get(self) -> float:
        """ Return a random value from gamma distribution.

        Examples:
        >>> for _ in range(100):
        ...     x = Gamma(1, 0, 0.627).get()
        ...     assert isinstance(x, float)
        ...     assert 0 <= x
        """
        a = self.a.get()
        loc = self.loc.get()
        scale = self.scale.get()
        rv = stats.gamma.rvs(
            a=a,
            loc=loc,
            scale=scale,
            size=None,
        )

        return float(rv)

    def get_many(self, n: int) -> List[float]:
        """ Return n random values from the gamma distribution.

        Examples:
        >>> for x in Gamma(1, 1, 100).get_many(100):
        ...     assert isinstance(x, float)
        ...     assert 1 <= x
        """
        a = self.a.get()
        loc = self.loc.get()
        scale = self.scale.get()
        rv = stats.gamma.rvs(
            a=a,
            loc=loc,
            scale=scale,
            size=n,
        )

        return rv.tolist()

    def seq(self, n: int) -> List[float]:
        """ Return n quantiles from the gamma distribution.

        Examples:
        >>> for x in Gamma(1, 5, 1).seq(100):
        ...     assert isinstance(x, float)
        ...     assert 5 <= x
        """

        a = self.a.get()
        loc = self.loc.get()
        scale = self.scale.get()

        quantiles = np.linspace(0.01, 0.99, num=n)
        return (stats.gamma
                .ppf(quantiles, a=a, loc=loc, scale=scale)
                .astype(float)
                .tolist())


class Normal(PDF):

    def __init__(
        self,
        loc: DistributionIF,
        scale: DistributionIF
    ) -> None:
        """ Take a random float from the normal distribution. """

        if isinstance(loc, float):
            self.loc: DistributionIF = FloatConst(loc)
        elif isinstance(loc, int):
            self.loc = FloatConst(float(loc))
        else:
            self.loc = loc

        if isinstance(scale, float):
            self.scale: DistributionIF = FloatConst(scale)
        elif isinstance(loc, int):
            self.scale = FloatConst(float(scale))
        else:
            self.scale = scale
        return

    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}"
                "({repr(self.loc)}, {repr(self.scale)})")

    def get(self) -> float:
        """ Return a random value from normal distribution.

        Examples:
        >>> for _ in range(100):
        ...     x = Normal(0, 0.627).get()
        ...     assert isinstance(x, float)
        """
        loc = self.loc.get()
        scale = self.scale.get()
        rv = stats.norm.rvs(
            loc=loc,
            scale=scale,
            size=None,
        )

        return float(rv)

    def get_many(self, n: int) -> List[float]:
        """ Return n random values from the normal distribution.

        Examples:
        >>> for x in Normal(1, 100).get_many(100):
        ...     assert isinstance(x, float)
        """
        loc = self.loc.get()
        scale = self.scale.get()
        rv = stats.norm.rvs(
            loc=loc,
            scale=scale,
            size=n,
        )

        return rv.tolist()

    def seq(self, n: int) -> List[float]:
        """ Return n quantiles from the normal distribution.

        Examples:
        >>> for x in Normal(5, 1).seq(100):
        ...     assert isinstance(x, float)
        """

        loc = self.loc.get()
        scale = self.scale.get()

        quantiles = np.linspace(0.01, 0.99, num=n)
        return (stats.norm
                .ppf(quantiles, loc=loc, scale=scale)
                .astype(float)
                .tolist())


class Range(PMF):

    def __init__(self, min: Distribution[int], max: Distribution[int]) -> None:
        """ Take a random int between min and max. """

        if isinstance(min, int):
            self.min: Distribution[int] = IntConst(min)
        else:
            self.min = min

        if isinstance(max, int):
            self.max: Distribution[int] = IntConst(max)
        else:
            self.max = max
        return

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({repr(self.min)}, {repr(self.max)})"

    def get(self) -> int:
        """ Return a random integer between min and max.

        Examples:
        >>> for _ in range(100):
        ...     x = Range(0, 5).get()
        ...     assert isinstance(x, int)
        ...     assert 0 <= x <= 5
        """

        min_ = self.min.get()
        max_ = self.max.get() + 1
        return stats.randint.rvs(min_, max_)

    def get_many(self, n: int) -> List[int]:
        """ Return n random values between min and max.

        Examples:
        >>> for x in Range(0, 10).get_many(100):
        ...     assert isinstance(x, int)
        ...     assert 0 <= x <= 10
        """
        min_ = self.min.get()
        max_ = self.max.get() + 1
        rv = stats.randint.rvs(min_, max_, size=n)
        return rv.tolist()

    def seq(self, n: int) -> List[int]:
        """ Return up to n unique quantiles between min and max.

        Examples:
        >>> for x in Range(0, 5).seq(100):
        ...     assert isinstance(x, int)
        ...     assert 0 <= x <= 5
        """

        min_ = self.min.get()
        max_ = self.max.get() + 1
        li = list(range(min_, max_))
        return self._find_boundaries(li, n)

    @staticmethod
    def _inner_boundaries(li: List[int], p: int) -> List[int]:
        boundaries = list()
        li = list(li)

        while p > 1 and len(li) > 0:
            i = len(li) // p
            boundaries.append(li[i])
            p -= 1
            li = li[i + 1:]
        return boundaries

    @classmethod
    def _find_boundaries(cls, li: List[int], n: int) -> List[int]:
        if len(li) <= n:
            return li
        elif n == 1:
            return [li[len(li) // 2]]

        head = li[0]
        tail = li[-1]

        n -= 1
        if n < 2:
            return [head, tail]
        else:
            return [head] + cls._inner_boundaries(li[1: -1], n) + [tail]


class Binomial(PMF):

    def __init__(
        self,
        n: Distribution[int],
        p: Distribution[float],
        loc: Distribution[int] = IntConst(0),
    ) -> None:
        """ Take a random integer from a binomial distribution. """

        if isinstance(n, int):
            self.n: Distribution[int] = IntConst(n)
        else:
            self.n = n

        if isinstance(loc, int):
            self.loc: Distribution[int] = IntConst(loc)
        else:
            self.loc = loc

        if isinstance(p, float):
            self.p: DistributionIF = FloatConst(p)
        else:
            self.p = p
        return

    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}"
                f"({repr(self.n)}, {repr(self.p)}, {repr(self.loc)})")

    def get(self) -> int:
        """ Return a random integer from the binomial distribution.

        Examples:
        >>> for _ in range(100):
        ...     x = Binomial(3, 0.5, 0).get()
        ...     assert isinstance(x, int)
        """

        n = self.n.get()
        loc = self.loc.get()
        p = self.p.get()
        return stats.binom.rvs(n, p, loc=loc, size=None)

    def get_many(self, n: int) -> List[int]:
        """ Return n random integers from a binomial distribution.

        Examples:
        >>> for x in Binomial(3, 0.5, 1).get_many(100):
        ...     assert isinstance(x, int)
        ...     assert 1 <= x
        """

        n_ = self.n.get()
        loc = self.loc.get()
        p = self.p.get()
        rv = stats.binom.rvs(n_, p, loc=loc, size=n)
        return rv.tolist()

    def seq(self, n: int) -> List[int]:
        """ Return up to n unique quantiles from a binomial distribution.

        Examples:
        >>> for x in Binomial(3, 0.2, 1).seq(10):
        ...     assert isinstance(x, int), "not int"
        ...     assert 1 <= x <= 4, "domain wrong"
        """

        n_ = self.n.get()
        loc = self.loc.get()
        p = self.p.get()
        quantiles = np.linspace(0.01, 0.99, n)
        return (np.unique(stats.binom.ppf(quantiles, n_, p, loc=loc))
                .astype(int)
                .tolist())


class NegBinomial(PMF):

    def __init__(
        self,
        n: Distribution[int],
        p: Distribution[float],
        loc: Distribution[int] = IntConst(0),
    ) -> None:
        """ Take a random int from the negative binomial distribution. """

        if isinstance(n, int):
            self.n: Distribution[int] = IntConst(n)
        elif isinstance(n, float):
            raise ValueError(n)
        else:
            self.n = n

        if isinstance(loc, int):
            self.loc: Distribution[int] = IntConst(loc)
        else:
            self.loc = loc

        if isinstance(p, float):
            self.p: Distribution[float] = FloatConst(p)
        else:
            self.p = p
        return

    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}"
                f"({repr(self.n)}, {repr(self.p)}, {repr(self.loc)})")

    def get(self) -> int:
        """ Return a random integer from the negative binomial distribution.

        Examples:
        >>> for _ in range(100):
        ...     x = NegBinomial(3, 0.5, 0).get()
        ...     assert isinstance(x, int)
        ...     assert x >= 0
        """

        n = self.n.get()
        loc = self.loc.get()
        p = self.p.get()
        return stats.nbinom.rvs(n, p, loc=loc, size=None)

    def get_many(self, n: int) -> List[int]:
        """ Return n random values from the negative binomial distribution.

        Examples:
        >>> for x in NegBinomial(3, 0.5, 1).get_many(100):
        ...     assert isinstance(x, int)
        ...     assert 1 <= x
        """

        n_ = self.n.get()
        loc = self.loc.get()
        p = self.p.get()
        rv = stats.nbinom.rvs(n_, p, loc=loc, size=n)
        return rv.tolist()

    def seq(self, n: int) -> List[int]:
        """ Return n random values from the negative binomial distribution.

        Examples:
        >>> for x in NegBinomial(3, 0.2, 1).seq(10):
        ...     assert isinstance(x, int), "not int"
        ...     assert 1 <= x
        """

        n_ = self.n.get()
        loc = self.loc.get()
        p = self.p.get()
        quantiles = np.linspace(0.01, 0.99, n)
        return (np.unique(stats.nbinom.ppf(quantiles, n_, p, loc=loc))
                .astype(int)
                .tolist())


class HyperGeometric(PMF):

    def __init__(
        self,
        n_success: Distribution[int],
        n: Distribution[int],
        population: Distribution[int],
        loc: Distribution[int] = IntConst(0),
    ) -> None:
        """ Take a random int from the hypergeometric distribution.
        """

        if isinstance(n_success, int):
            self.n_success: Distribution[int] = IntConst(n_success)
        elif isinstance(n_success, float):
            raise ValueError(n_success)
        else:
            self.n_success = n_success

        if isinstance(n, int):
            self.n: Distribution[int] = IntConst(n)
        elif isinstance(n, float):
            raise ValueError(n)
        else:
            self.n = n

        if isinstance(population, int):
            self.population: Distribution[int] = IntConst(population)
        elif isinstance(population, float):
            raise ValueError(population)
        else:
            self.population = population

        if isinstance(loc, int):
            self.loc: Distribution[int] = IntConst(loc)
        else:
            self.loc = loc

        return

    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}"
                f"({repr(self.n_success)}, {repr(self.n)}, "
                f"{repr(self.population)}, {repr(self.loc)})")

    def get(self) -> int:
        """ Return a random integer from the hypergeometric distribution.

        Examples:
        >>> for _ in range(100):
        ...     x = HyperGeometric(7, 12, 20, 0).get()
        ...     assert isinstance(x, int)
        ...     assert x >= 0
        """

        n_success = self.n_success.get()
        n = self.n.get()
        population = self.population.get()
        loc = self.loc.get()
        return stats.hypergeom.rvs(
            M=population,
            n=n_success,
            N=n,
            loc=loc,
            size=None
        )

    def get_many(self, n: int) -> List[int]:
        """ Return n random values from the hypergeometric distribution.

        Examples:
        >>> for x in HyperGeometric(7, 12, 20, 1).get_many(100):
        ...     assert isinstance(x, int)
        ...     assert 1 <= x
        """

        n_success = self.n_success.get()
        n_ = self.n.get()
        population = self.population.get()
        loc = self.loc.get()
        rv = stats.hypergeom.rvs(
            M=population,
            n=n_success,
            N=n_,
            loc=loc,
            size=n
        )
        return rv.tolist()

    def seq(self, n: int) -> List[int]:
        """ Return n quantiles from the hypergeometric distribution.

        Examples:
        >>> for x in HyperGeometric(7, 12, 20, 1).seq(10):
        ...     assert isinstance(x, int), "not int"
        ...     assert 1 <= x
        """

        n_success = self.n_success.get()
        n_ = self.n.get()
        population = self.population.get()
        loc = self.loc.get()

        quantiles = np.linspace(0.01, 0.99, n)
        rvs = stats.hypergeom.ppf(
            quantiles,
            population,
            n_success,
            n_,
            loc=loc
        )
        return np.unique(rvs).astype(int).tolist()


class Choose(Generic[T]):

    def __init__(self, choices: Sequence[Distribution[T]]) -> None:
        """ Choose an element from a list of possible choices. """

        self.choices: List[Distribution[T]] = list(choices)
        return

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({repr(self.choices)})"

    def get(self) -> T:
        """ Returns a random choice. """

        import random
        choice = random.choice(self.choices)
        return choice.get()

    def get_many(self, n: int) -> List[T]:
        """ Returns n random choices. """

        import random
        choices = random.choices(self.choices, k=n)
        return [c.get() for c in choices]


class ChooseI(Choose[int], PMF):
    """ Concrete type for Choose[int] """
    pass


class ChooseF(Choose[Union[float, int]], PDF):
    """ Concrete type for Choose[float] """
    pass


class ChooseS(Choose[str], Categorical):
    """ Concrete type for Choose[str] """
    pass


class Trunc(Distribution[V]):

    def __init__(self, min: V, max: V, dist: Distribution[V]) -> None:
        """ A wrapper to prevent a distribution from selecting an invalid
        number.

        Mostly for enforcing X >= 0 in normal distributions.
        Selecting one of these values should be a rare event.
        """

        self.min: V = min
        self.max: V = max
        self.dist: Distribution[V] = dist
        return

    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}"
                "({self.min}, {self.max}, {repr(self.dist)})")

    def get(self) -> V:
        """ Get a sample from the distribution within the bounds.

        Raises `DistDomainError` if haven't found suitable value within
        a number of attempts.

        Examples:
        >>> dist = Uniform(0, 2)
        >>> tdist = Trunc(1, 1.5, dist)
        >>> for _ in range(100):
        ...     x = tdist.get()
        ...     assert 1 <= x <= 1.5
        """

        max_attempts = 1000
        for _ in range(max_attempts):
            attempt = self.dist.get()
            if self.min <= attempt <= self.max:
                return attempt

        raise DistDomainError(
            f"Could not get value between {self.min} and {self.max} within "
            f"{max_attempts} from distribution {repr(self.dist)}. "
            f"Please ensure that the distribution specification is reasonable "
            "given the truncation bounds."
        )

    def get_many(self, n: int) -> List[V]:
        """ Get n samples from the distribution within the bounds.

        Raises `DistDomainError` if haven't found suitable value within
        a number of attempts.

        Examples:
        >>> dist = Uniform(0, 2)
        >>> tdist = Trunc(1, 1.5, dist)
        >>> for x in tdist.get_many(100):
        ...     assert 1 <= x <= 1.5
        """

        max_attempts = 1000
        for _ in range(max_attempts):
            attempt = self.dist.get_many(n * max_attempts)
            filtered_attempt = [
                a
                for a
                in attempt
                if self.min <= a <= self.max
            ]

            if len(filtered_attempt) >= n:
                return filtered_attempt[:n]

        raise DistDomainError(
            f"Could not get n values between {self.min} and {self.max} within "
            f"{max_attempts} from distribution {repr(self.dist)}. "
            "Please ensure that the distribution specification is reasonable "
            "given the truncation bounds."
        )

    def seq(self, n: int) -> List[V]:
        """ Drops any elements outside the bounds.

        Examples:
        >>> tdist = Trunc(1, 2, Uniform(0, 2))
        >>> tdist.seq(11)
        [1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
        """
        # Seq should always be monotonically increasing
        seq = self.dist.seq(n)
        fseq = [s for s in seq if self.min <= s <= self.max]

        if len(fseq) == 0:
            raise DistDomainError(
                f"Could not get any values from quantiles between {self.min} "
                f"and {self.max} from distribution {repr(self.dist)}. "
                f"Values were {seq}. "
                "Please ensure that the distribution specification is "
                "reasonable given the truncation bounds."
            )

        return fseq


# Sequence types


class AOList(Generic[T]):

    def __init__(self, elems: Sequence[Distribution[T]]) -> None:
        self.elems: List[Distribution[T]] = list(elems)
        return

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({repr(self.elems)})"

    def get(self) -> List[T]:
        return [e.get() for e in self.elems]

    def __iter__(self) -> Iterable[T]:
        return map(lambda x: x.get(), iter(self.elems))

    def __len__(self) -> int:
        return len(self.elems)


class Sample(AOList[T]):

    def __init__(self, n: int, dist: Distribution[T]) -> None:
        self.n = n
        self.dist: Distribution[T] = dist
        return

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.n}, {repr(self.dist)})"

    def get(self) -> List[T]:
        return self.dist.get_many(self.n)


class Monotonic(AOList[T]):

    def __init__(self, n: int, dist: Distribution[T]) -> None:
        """ Get a monotonically increasing list from a distribution.

        Example:
        >>> m = Monotonic(10, Range(0, 15))
        >>> m
        Monotonic(10, Range(IntConst(0), IntConst(15)))
        >>> vals = m.get()
        >>> prev = vals[0]
        >>> for v in vals[1:]:
        ...     assert v > prev
        ...     prev = v
        >>> m = Monotonic(20, Range(1, 10))
        >>> vals = m.get()
        >>> assert len(vals) == 10
        >>> m = Monotonic(20, Uniform(1, 10))
        >>> vals = m.get()
        >>> assert len(vals) == 20
        """

        self.n = n
        self.dist: Distribution[T] = dist
        return

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.n}, {repr(self.dist)})"

    def get(self) -> List[T]:
        seen: Set[T] = set()

        i = 0
        while (len(seen) < self.n) and (i < self.n + 1000):
            val = self.dist.get()
            if val not in seen:
                seen.add(val)

            i += 1

        return sorted(seen)


class Seq(AOList[T]):

    def __init__(self, n: int, dist: Distribution[T]) -> None:
        """ """
        self.n: int = n
        self.dist: Distribution[T] = dist
        return

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.n}, {repr(self.dist)})"

    def get(self) -> List[T]:
        """ Returns a partition over the quantiles of the distribution. """

        return self.dist.seq(self.n)

# seq 5 uniform 0.5 0.9
# 0.5 0.6 0.7 0.8 0.9
