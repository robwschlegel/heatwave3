"""Thread-count control (port of R/zzz.R)."""

from . import _core


def getHW3threads():
    """Number of threads heatwave3 uses by default.

    Defaults to half the available cores, or the ``HEATWAVE3_NUM_THREADS``
    (or ``R_HEATWAVE3_NUM_THREADS``) environment variable when set.
    """
    return _core.hw3_get_threads()


def setHW3threads(threads=0):
    """Set the default number of threads; ``0`` restores the default.

    Returns
    -------
    int
        The number of threads now in effect.
    """
    return _core.hw3_set_threads(int(threads))
