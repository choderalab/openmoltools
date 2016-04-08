#!/usr/bin/env python

"""Miscellaneous wrapper functions to Schrodinger's computational chemistry tools."""

import os
import sys
import csv
import shutil
import logging
import subprocess
from functools import wraps

import mdtraj

from openmoltools import utils

logger = logging.getLogger(__name__)


def run_and_log_error(command):
    """Run the process specified by the command and log eventual errors.

    Parameters
    ----------
    command : str
        The command to be run.

    Returns
    -------
    output : str
        The output of the process.

    Raises
    ------
    subprocess.CalledProcessError
        In case the commands fails.

    """
    try:
        output = subprocess.check_output(command)
    except subprocess.CalledProcessError as e:
        logger.error(e.output)
        logger.error(str(e))
        raise e
    return output.decode()


def is_schrodinger_suite_installed():
    """Check that Schrodinger's suite is installed.

    Currently only checks whether the environmental variable SCHRODINGER
    is defined. This should contain the path to its main installation folder.

    Returns
    -------
    bool
        True if the Schrodinger's suite is found, False otherwise.

    """
    try:
        os.environ['SCHRODINGER']
    except KeyError:
        return False
    return True


def need_schrodinger(func):
    @wraps(func)
    def _need_schrodinger(*args, **kwargs):
        """Decorator that checks if the Schrodinger's suite is installed."""
        if not is_schrodinger_suite_installed():
            err_msg = "Cannot locate Schrodinger's suite!"
            logger.error(err_msg)
            raise RuntimeError(err_msg)
        return func(*args, **kwargs)
    return _need_schrodinger


@need_schrodinger
def run_proplister(input_file_path):
    """Run proplister utility on a file and return its properties.

    Parameters
    ----------
    input_file_path: str
        The path to the file describing the molecule with its properties.

    Returns
    -------
    properties: list of dict
        A list containing a dictionary for each molecule in the input file
        representing their properties. Each dictionary is in the format
        property_name -> property_value.

    """
    proplister_path = os.path.join(os.environ['SCHRODINGER'], 'utilities', 'proplister')

    # Normalize path
    input_file_path = os.path.abspath(input_file_path)

    # Run proplister, we need the list in case there are spaces in paths
    cmd = [proplister_path, '-a', '-c', input_file_path]
    output = run_and_log_error(cmd)

    output = output.replace('\_', '_')  # Parse '\_' characters in names

    # The output is a cvs file. The first line are the property names and then each row
    # contains the values for each molecule. We use the csv module to avoid splitting
    # strings that contain commas (e.g. "2,2-dimethylpropane").
    properties = []
    csv_reader = csv.reader(output.split('\n'))
    names = next(csv_reader)
    for values in csv_reader:
        if len(values) == 0:
            continue  # proplister prints a final empty line

        # Convert raw strings into literals (e.g. convert '\\n' to '\n')
        if sys.version_info < (3, 0):  # Python 2
            converted_values = [v.decode('string_escape') for v in values]
        else:  # Python 3 doesn't have decode on strings
            converted_values = [bytes(v, "utf-8").decode("unicode_escape")
                                for v in values]

        properties.append(dict(zip(names, converted_values)))

    return properties


@need_schrodinger
def run_structconvert(input_file_path, output_file_path):
    """Run Schrodinger's structconvert command line utility to convert from one
    format to another.

    The input and output formats are inferred from the given files extensions.

    Parameters
    ----------
    input_file_path : str
        Path to the input file describing the molecule.
    output_file_path : str
        Path were the converted file will be saved.

    """
    formats_map = {'sdf': 'sd'}  # convert common extensions to format code

    # Locate structconvert executable
    structconvert_path = os.path.join(os.environ['SCHRODINGER'], 'utilities',
                                      'structconvert')

    # Normalize paths
    input_file_path = os.path.abspath(input_file_path)
    output_file_path = os.path.abspath(output_file_path)

    # Determine input and output format
    input_format = os.path.splitext(input_file_path)[1][1:]
    output_format = os.path.splitext(output_file_path)[1][1:]
    if input_format in formats_map:
        input_format = formats_map[input_format]
    if output_format in formats_map:
        output_format = formats_map[output_format]

    # Run structconvert, we need the list in case there are spaces in paths
    cmd = [structconvert_path, '-i' + input_format, input_file_path,
           '-o' + output_format, output_file_path]
    run_and_log_error(cmd)


def autoconvert_maestro(func):
    @wraps(func)
    def _autoconvert_maestro(input_file_path, output_file_path, *args, **kwargs):
        """Decorator that make a function support more than only Maestro files.

        Input and output formats are inferred from extensions. If the input file
        is not in Maestro format, this automatically uses the utility structconvert
        to create a temporary Maestro file. Similarly, if the output file path does
        not have a 'mae' extension, a temporary output file is created and converted
        at the end of the wrapped function execution.

        The decorated function must take as first two parameters the input and the
        output paths respectively.

        """
        is_input_mae = os.path.splitext(input_file_path)[1] == '.mae'
        is_output_mae = os.path.splitext(output_file_path)[1] == '.mae'

        # If they are both in Maestro format just call the function
        if is_output_mae and is_input_mae:
            return func(input_file_path, output_file_path, *args, **kwargs)

        # Otherwise we create a temporary directory to host temp files
        # First transform desired paths into absolute
        input_file_path = os.path.abspath(input_file_path)
        output_file_path = os.path.abspath(output_file_path)
        with mdtraj.utils.enter_temp_directory():
            # Convert input file if necessary
            if is_input_mae:
                func_input = input_file_path
            else:
                func_input = os.path.splitext(os.path.basename(input_file_path))[0] + '.mae'
                run_structconvert(input_file_path, func_input)

            # Determine if we need to convert output
            if is_output_mae:
                func_output = output_file_path
            else:
                func_output = os.path.splitext(os.path.basename(output_file_path))[0] + '.mae'

            # Execute function
            return_value = func(func_input, func_output, *args, **kwargs)

            # Delete temporary input
            if not is_input_mae:
                os.remove(func_input)

            # Convert temporary output
            if not is_output_mae:
                run_structconvert(func_output, output_file_path)
                os.remove(func_output)

            # Copy any other output file in the temporary folder
            output_dir = os.path.dirname(output_file_path)
            for file_name in os.listdir('.'):
                shutil.copy2(file_name, os.path.join(output_dir, file_name))

        return return_value
    return _autoconvert_maestro


@need_schrodinger
@autoconvert_maestro
def run_maesubset(input_file_path, output_file_path, range):
    """Run Schrodinger's maesubset command line utility to extract a range of
    structures from a file.

    Parameters
    ----------
    input_file_path : str
        Path to the input file with multiple structures.
    output_file_path : str
        Path to output file.
    range : int or list of ints
        The 0-based indices of the structures to extract from the input files.

    """

    # Locate maesubset executable
    maesubset_path = os.path.join(os.environ['SCHRODINGER'], 'utilities', 'maesubset')

    # Normalize paths
    input_file_path = os.path.abspath(input_file_path)
    output_file_path = os.path.abspath(output_file_path)

    # Determine molecules to extract
    try:  # if range is a list of ints
        range_str = [str(i + 1) for i in range]
    except TypeError:  # if range is an int
        range_str = [str(range + 1)]
    range_str = ','.join(range_str)

    # Run maesubset, we need the list in case there are spaces in paths
    cmd = [maesubset_path, '-n', range_str, input_file_path]
    output = run_and_log_error(cmd)

    # Save result
    with open(output_file_path, 'w') as f:
        f.write(output)


@need_schrodinger
@autoconvert_maestro
def run_epik(input_file_path, output_file_path, max_structures=32, ph=7.4,
             ph_tolerance=None, min_probability=None, tautomerize=True, extract_range=None):
    """Run Schrodinger's epik command line utility to enumerate protonation and
    tautomeric states.

    Parameters
    ----------
    input_file_path : str
        Path to input file describing the molecule.
    output_file_path : str
        Path to the output file created by epik.
    max_structures : int, optional
        Maximum number of generated structures (default is 32).
    ph : float, optional
        Target pH for generated states (default is 7.4).
    ph_tolerance : float, optional
        Equivalent of -pht option in Epik command (default is None).
    min_probability: float, optional
        Minimum probability for the generated states.
    tautomerize : bool, optional
        Whether or not tautomerize the input structure (default is True).
    extract_range : int or list of ints, optional
        If not None, the function uses the Schrodinger's utility maesubset to
        extract only a subset of the generated structures. This is the 0-based
        indices of the structures to extract from the input files.
    """

    # Locate epik executable
    epik_path = os.path.join(os.environ['SCHRODINGER'], 'epik')

    # Normalize paths as we'll run in a different working directory
    input_file_path = os.path.abspath(input_file_path)
    output_file_path = os.path.abspath(output_file_path)
    output_dir = os.path.dirname(output_file_path)

    # Preparing epik command arguments for format()
    epik_args = dict(ms=max_structures, ph=ph)
    epik_args['pht'] = '-pht {}'.format(ph_tolerance) if ph_tolerance else ''
    epik_args['nt'] = '' if tautomerize else '-nt'
    epik_args['p'] = '-p {}'.format(min_probability) if min_probability else ''

    # Determine if we need to convert input and/or output file
    if extract_range is None:
        epik_output = output_file_path
    else:
        epik_output = os.path.splitext(output_file_path)[0] + '-full.mae'

    # Epik command. We need list in case there's a space in the paths
    cmd = [epik_path, '-imae', input_file_path, '-omae', epik_output]
    cmd += '-ms {ms} -ph {ph} {pht} {nt} {p} -pKa_atom -WAIT -NO_JOBCONTROL'.format(
            **epik_args).split()

    # We run with output_dir as working directory to save there the log file
    with utils.temporary_cd(output_dir):
        run_and_log_error(cmd)

    # Check if we need to extract a range of structures
    if extract_range is not None:
        run_maesubset(epik_output, output_file_path, extract_range)
        os.remove(epik_output)
