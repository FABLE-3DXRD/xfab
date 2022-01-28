import logging

def get_module_level_logger(module):
    """Build and set up formatting for a module level logger.

    This is meant to be called in each xfab module as:

        logger = get_module_level_logger(__name__)

    providing a module level logger to use for all local logging.

    The reason for module level loggers is to allow the end user to see what modules produced the logs
    and give controll over the level of logging. For instance,

        logging.getLogger('xfab.structure').setLevel(CRITICAL)

    will only display logging messages which are errors (or above in severity). Likewise following will shut down
    logging from xfab.structure completely:

        logging.getLogger('xfab.structure').disabled = True

    Some more logger docs can be found here: https://docs.python.org/3/howto/logging.html#configuring-logging

    Args:
        module (:obj:`str`): The module name and parents. e.g 'xfab.structure'.

    Returns:
        (:obj:`logging.logger`) A module level logger.

    """
    # create logger on lowest level=NOTSET
    logger = logging.getLogger(module)
    logger.setLevel(logging.NOTSET)

    # create console handler and set level
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.NOTSET)

    # Some formatting printing what module is logging and the severity etc
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # add the formatter to the console handler
    console_handler.setFormatter(formatter)

    # add the console handler to logger
    logger.addHandler(console_handler)

    return logger
