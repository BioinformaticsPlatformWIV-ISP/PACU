from pathlib import Path


def path_to_absolute(path: str) -> Path:
    """
    Takes a relative path and returns the absolute path.
    :param path: Relative or absolute path
    :return: The resolved absolute Path object
    """
    return Path(path).expanduser().resolve()
