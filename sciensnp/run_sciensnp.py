#!/usr/bin/env python
from sciensnp import initialize_logging, ScienSNP


def main() -> None:
    """
    Runs the main workflow.
    :return: None
    """
    initialize_logging()
    main_ = ScienSNP()
    main_.run()


if __name__ == '__main__':
    main()
