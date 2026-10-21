"""Centralized logging configuration for PREFACE."""

from __future__ import annotations

import logging
import sys
from enum import Enum
from typing import TextIO

LOGGER_NAME = "preface"


class LogLevel(str, Enum):
    """Supported CLI verbosity levels."""

    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"
    CRITICAL = "CRITICAL"


class _StderrHandler(logging.StreamHandler[TextIO]):
    """StreamHandler that dynamically resolves the current sys.stderr stream."""

    @property
    def stream(self) -> TextIO:
        return sys.stderr

    @stream.setter
    def stream(self, value: TextIO | None) -> None:
        pass


def configure_logging(level: str | LogLevel = LogLevel.INFO) -> logging.Logger:
    """Configure the package logger to emit formatted messages to stderr."""
    level_name = (
        level.value if isinstance(level, LogLevel) else str(level).strip().upper()
    )
    if level_name not in LogLevel.__members__:
        valid = ", ".join(LogLevel.__members__.keys())
        raise ValueError(f"Unsupported log level '{level}'. Choose from: {valid}.")

    numeric_level = getattr(logging, level_name)
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(numeric_level)
    logger.propagate = False

    fmt = (
        "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"
        if numeric_level == logging.DEBUG
        else "%(levelname)s: %(message)s"
    )
    formatter = logging.Formatter(fmt=fmt, datefmt="%Y-%m-%d %H:%M:%S")

    if not logger.handlers:
        handler = _StderrHandler()
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    else:
        for existing in logger.handlers:
            existing.setFormatter(formatter)

    return logger


def get_logger(name: str = LOGGER_NAME) -> logging.Logger:
    """Return a child logger under the preface namespace."""
    if name == LOGGER_NAME or name.startswith(f"{LOGGER_NAME}."):
        return logging.getLogger(name)
    return logging.getLogger(f"{LOGGER_NAME}.{name}")
