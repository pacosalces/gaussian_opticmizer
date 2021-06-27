import os, sys
import yaml


def loader(part_no):
    """Loads params from part_no.yml file as dict

    Args:
        part_no (str): filename of requested optic
    """
    with open(part_no + ".yaml", "r") as lens_datasheet:
        return yaml.safe_load(lens_datasheet)


if __name__ == "__main__":
    loader("./LA4380")
