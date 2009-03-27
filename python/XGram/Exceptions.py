class Error(Exception):
    """Base class for exceptions in this module."""
    pass

    def __str__(self):
        return str(self.message)

class ParsingError(Error):
    """Exception raised for errors while parsing

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message, line):
        self.message = message + " at line " + line

class UsageError(Error):
    """Exception raised for errors while starting

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class ProgrammingError(Error):
    """Exception raised for errors that result from programming errors.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
