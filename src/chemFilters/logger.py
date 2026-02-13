# -*- coding: utf-8 -*-

import sys

from loguru import logger

COLORFUL_FORMAT = (
    "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}"
    "</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan>"
    " - <level>{message}</level>"
)


def setup_logger(level="INFO"):
    """Add a colorful log handler filtered to chemFilters messages only.

    Does not remove existing handlers. Returns the handler ID so it can be
    removed later with ``logger.remove(handler_id)``.
    """
    handler_id = logger.add(
        sys.stderr,
        format=COLORFUL_FORMAT,
        level=level,
        filter="chemFilters",
        colorize=True,
    )
    return handler_id
