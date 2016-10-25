"""
Command line interface for mdgenerate.
"""

from .mdgenerate import process

import argparse
import os


def run(args=None):
    """
    Run the mdgenerate command line interface.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        'yaml',
        help='The YAML file to parse.',
        nargs='?',
        default='meta.yaml'
    )
    grompp_parser = parser.add_mutually_exclusive_group(required=False)
    grompp_parser.add_argument(
        '--grompp',
        help='If grompp should be run',
        dest='grompp', action='store_true'
    )
    grompp_parser.add_argument(
        '--no-grompp',
        dest='grompp', action='store_false'
    )
    parser.set_defaults(grompp=True)
    args = parser.parse_args()

    yaml = os.path.abspath(args.yaml)

    if os.path.exists(yaml):
        process(yaml, run_grompp=args.grompp)
    else:
        raise FileNotFoundError('File not found: {}'.format(yaml))
