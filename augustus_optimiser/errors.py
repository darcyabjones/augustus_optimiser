from enum import IntEnum
from typing import ClassVar


class ECode(IntEnum):
    """ Error codes for command termination
    Based on linux '/usr/include/sysexits.h'
    Except for usage errors which are 2, in keeping with argparse.
    """

    OK = 0
    ERROR = 1  # Generic error, avoid using if possible
    USAGE = 2  # command line usage error
    DATAERR = 65  # data format error
    NOINPUT = 66  # cannot open input
    NOUSER = 67  # addressee unknown
    NOHOST = 68  # host name unknown
    UNAVAILABLE = 69  # service unavailable
    SOFTWARE = 70  # internal software error
    OSERR = 71  # system error (e.g., can't fork)
    OSFILE = 72  # critical OS file missing
    CANTCREAT = 73  # can't create (user) output file
    IOERR = 74  # input/output error
    TEMPFAIL = 75  # temp failure; user is invited to retry
    PROTOCOL = 76  # remote error in protocol
    NOPERM = 77  # permission denied
    CONFIG = 78  # configuration error
    SIGINT = 130  # Ctrl-c


class AOError(Exception):

    ecode: ClassVar[ECode] = ECode.ERROR

    def __init__(self, msg: str) -> None:
        self.msg = msg
        return

    def __str__(self) -> str:
        return f"{self.msg}"


class DistParseError(AOError):

    ecode: ClassVar[ECode] = ECode.DATAERR


class DistValueError(AOError):

    ecode: ClassVar[ECode] = ECode.SOFTWARE


class MaxAttemptsError(AOError):

    ecode: ClassVar[ECode] = ECode.SOFTWARE


class DistDomainError(AOError):

    ecode: ClassVar[ECode] = ECode.SOFTWARE


class ConfigParseError(AOError):

    ecode: ClassVar[ECode] = ECode.DATAERR
