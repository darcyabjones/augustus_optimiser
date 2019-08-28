import sys
from ao.cli import cli_not_subparser, main_error_handler
from ao.extrinsic import run_extrinsic, cli_extrinsic, EXTRINSIC_DESC


def main():
    args = cli_not_subparser(
        prog="ao-extrinsic",
        args=sys.argv[1:],
        desc=EXTRINSIC_DESC,
        subparser=cli_extrinsic
    )
    main_error_handler(args, run_extrinsic)
    return


if __name__ == "__main__":
    main()
