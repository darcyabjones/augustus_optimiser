import sys
import argparse


from ao.errors import AOError, ECode
from ao.extrinsic import run_extrinsic, cli_extrinsic, EXTRINSIC_DESC
from ao.extract import run_extract, cli_extract, EXTRACT_DESC


def cli(prog, args):

    parser = argparse.ArgumentParser(
        prog=prog,
        description=""
    )

    subparsers = parser.add_subparsers(dest='subparser_name')

    extrinsic_subparser = subparsers.add_parser(
        "extrinsic",
        help=EXTRINSIC_DESC
    )

    cli_extrinsic(extrinsic_subparser)

    extract_subparser = subparsers.add_parser(
        "extract",
        help=EXTRACT_DESC
    )

    cli_extract(extract_subparser)

    parsed = parser.parse_args(args)

    if parsed.subparser_name is None:
        parser.print_help()
        sys.exit(0)

    return parsed


def cli_not_subparser(prog, args, desc, subparser):

    parser = argparse.ArgumentParser(
        prog=prog,
        description=desc
    )

    subparser(parser)
    parsed = parser.parse_args(args)

    return parsed


def main_error_handler(args, runner):
    try:
        runner(args)
    except AOError as e:
        print(f"Error: {e.msg}")
        sys.exit(e.ecode)
    except BrokenPipeError:
        # Pipes get closed and that's normal
        sys.exit(0)
    except KeyboardInterrupt:
        print("Received keyboard interrupt. Exiting.", file=sys.stderr)
        sys.exit(ECode.SIGINT)
    except EnvironmentError as e:
        print((
            "Encountered a system error.\n"
            "We can't control these, and they're usually related to your OS.\n"
            "Try running again."
        ), file=sys.stderr)
        raise e
    except Exception as e:
        print((
            "I'm so sorry, but we've encountered an unexpected error.\n"
            "This shouldn't happen, so please file a bug report with the "
            "authors.\nWe will be extremely grateful!\n\n"
        ), file=sys.stderr)
        raise e
    return


def main():
    args = cli(prog="ao", args=sys.argv[1:])
    if args.subparser_name == "extrinsic":
        main_error_handler(args, run_extrinsic)
    if args.subparser_name == "extract":
        main_error_handler(args, run_extract)
    else:
        raise ValueError("I shouldn't reach this point ever")
    return
