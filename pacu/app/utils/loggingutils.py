import logging


logger = logging.getLogger('pacu')


def initialize_logging() -> None:
    """
    Initializes the logging.
    :return: None
    """
    formatter = logging.Formatter('%(asctime)s - %(module)15s - %(levelname)7s - %(message)s')

    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    console_handler.setLevel(logging.DEBUG)
    console_handler.name = 'console'
    logger.addHandler(console_handler)

    # File handler (camel.log file)
    file_handler = logging.FileHandler('camel.log')
    file_handler.setFormatter(formatter)
    file_handler.setLevel(logging.DEBUG)
    logger.addHandler(file_handler)

    # General logging level
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
