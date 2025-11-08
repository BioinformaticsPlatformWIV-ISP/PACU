#!/usr/bin/env python
from pacu import PACU
from pacu.app.utils.loggingutils import initialize_logging


def main() -> None:
    """
    Runs the main workflow.
    :return: None
    """
    initialize_logging()
    main_ = PACU()
    main_.run()


if __name__ == '__main__':
    main()
