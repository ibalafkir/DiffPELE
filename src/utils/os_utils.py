import os
import subprocess as sp
from contextlib import contextmanager
import pandas as pd


def execute_command(cmd: str) -> str:
    """
    Execute a command in the shell and return the output
    :param cmd: The command to execute as a string.
    :return: The output of the command as a string.
    """
    # Split the command string into a list of arguments
    cmd_list = cmd.split()
    result = sp.run(cmd_list, capture_output=True, text=True)
    if result.returncode != 0:
        raise Exception(f"Command failed with error: {result.stderr}")
    return result.stdout

def display_full_dataframe():
    """
    Configures pandas to display complete DataFrames without truncation.
    """
    pd.set_option('display.max_rows', None)  # Show all rows
    pd.set_option('display.max_columns', None)  # Show all columns
    pd.set_option('display.width', None)  # Adjust width to avoid wrapping
    pd.set_option('display.max_colwidth', None)  # Show full column contents