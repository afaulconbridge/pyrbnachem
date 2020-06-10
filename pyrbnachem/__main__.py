import argparse
import logging
import logging.config

import random
import itertools
from .rbnachem import RBNAchem
from .rbnmol import RBNMol
import pyachem


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--log-config", type=str, default=None, help="log configuration file"
    )
    parser.add_argument(
        "--log-critical", action="store_const", const=logging.CRITICAL, dest="log_level"
    )
    parser.add_argument(
        "--log-error", action="store_const", const=logging.ERROR, dest="log_level"
    )
    parser.add_argument(
        "--log-warning", action="store_const", const=logging.WARNING, dest="log_level"
    )
    parser.add_argument(
        "--log-info",
        action="store_const",
        const=logging.INFO,
        dest="log_level",
        default=logging.INFO,
    )
    parser.add_argument(
        "--log-debug", action="store_const", const=logging.DEBUG, dest="log_level"
    )
    args = parser.parse_args()

    logging.basicConfig(level=args.log_level)
    if args.log_config:
        logging.config.fileConfig(args.log_config)

    # TODO update
    logger = logging.getLogger("TODO")

    if args.log_config:
        logger.debug("loaded log config from {}".format(args.log_config))

    # TODO finish

    achem = RBNAchem()
    contents = []
    for seed in range(10):
        mol = RBNMol.from_random(random.Random(seed), n=5, k=2)
        for _ in range(10):
            contents.append(mol)
    vessel = pyachem.VesselZip(achem, contents, random.Random(42))
    for reaction in itertools.islice(vessel, 10):
        print(reaction)


if __name__ == "__main__":
    main()
