class BMDSException(Exception):
    """An exception occurred in the BMDS module."""

    pass


class ConfigurationException(BMDSException):
    """The configuration of the BMDS module is invalid."""

    pass
