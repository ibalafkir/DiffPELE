import logging


class LoggerFactory:
    """
    A class to initialize the logger object
    Source: https://medium.com/geekculture/create-a-reusable-logger-factory-for-python-projects-419ad408665d
    """
    _LOG = None

    @staticmethod
    def _create_logger(log_file: str,
                       log_level: str):
        """
        A private method that interacts with the python
        logging module
        :param log_file:
            The name of the file from which the logger is called
        :param log_level:
            The level of logging. It can be INFO, ERROR, DEBUG
        """
        # set the logging format
        #log_format = "%(asctime)s:%(levelname)s:%(message)s"
        log_format = '%(module)s - %(levelname)s: %(message)s'

        # Initialize the class variable with logger object
        LoggerFactory._LOG = logging.getLogger(log_file)
        logging.basicConfig(level=logging.INFO, format=log_format,
                            datefmt="%Y-%m-%d %H:%M:%S")

        # set the logging level based on the user selection
        if log_level.upper() == "INFO":
            LoggerFactory._LOG.setLevel(logging.INFO)
        elif log_level.upper() == "ERROR":
            LoggerFactory._LOG.setLevel(logging.ERROR)
        elif log_level.upper() == "DEBUG":
            LoggerFactory._LOG.setLevel(logging.DEBUG)
        else:
            raise ValueError("Invalid logging level")
        return LoggerFactory._LOG

    @staticmethod
    def get_logger(log_file: str,
                   log_level: str):
        """
        A static method called by other modules to initialize logger in
        their own module
        """
        logger = LoggerFactory._create_logger(log_file, log_level)

        return logger