"""
Contains HotspotFrameworkerror class
"""


class HotspotFrameworkerror(Exception):
    """Base class for exceptions in HotspotFramework"""
    pass


class UserInputError(HotspotFrameworkerror):
    """
    Exception raised for errors in the input
    """
    def __init__(self, message):
        """
        Args:
            message (str): error found

        Returns:
            None
        """
        self.message = message
