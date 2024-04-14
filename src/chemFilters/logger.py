# -*- coding: utf-8 -*-
import sys

from loguru import logger


def setup_logger(level="INFO"):
    colorful_format = (
        "<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}"
        "</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan>"
        " - <level>{message}</level>"
    )
    logger.remove()
    logger.add(sys.stderr, format=colorful_format, level=level, colorize=True)


setup_logger()
