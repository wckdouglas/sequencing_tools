class SeqUtilsError(ValueError):
    """
    Exception type for the package sequencing tools

    Example::

        a = -1
        if a < 0:
           raise SeqUtilsError("a should be > 0")
    """
    pass
