"""
Used to store global variables that are used throughout the package.
Vars that can be set by the user are listed in the beginning of the file.

How to import and use:
import src.global_vars as gv
gv.variable_name

"""
import os


# User-settable variables

LOGGER_LEVEL = 'INFO'  # Set the logging level. Options: 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'


# Constants

AAs = 'ACDEFGHIKLMNPQRSTVWY'
AAs_3 = {
    'A': 'ALA',
    'C': 'CYS',
    'D': 'ASP',
    'E': 'GLU',
    'F': 'PHE',
    'G': 'GLY',
    'H': 'HIS',
    'I': 'ILE',
    'K': 'LYS',
    'L': 'LEU',
    'M': 'MET',
    'N': 'ASN',
    'P': 'PRO',
    'Q': 'GLN',
    'R': 'ARG',
    'S': 'SER',
    'T': 'THR',
    'V': 'VAL',
    'W': 'TRP',
    'Y': 'TYR'
}