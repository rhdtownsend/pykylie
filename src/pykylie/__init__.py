# -*- coding: utf-8 -*-

__version__ = "0.1+dev"

"""Python implementation of KYLIE.

PyKYLIE provides a set of routines that allow a simple Python
re-implementation of the KYLIE pulsating-star spectral synthesis code.

"""

import numpy as np
import re
import astropy.table as at
import astropy.constants as ac
import scipy.io as si
import pymsg as pm

#
    
def read_bruce_model(file_name):

    r"""Read BRUCE model data from a file.

    Args:
        file_name (str): Name of file.

    Returns:
        astropy.table.Table: Model data. Array values are stored in table
        columns, and scalar values in the table `meta` attribute.
    """

    with si.FortranFile(file_name, 'r') as f:
    
        # Read number of points and time
    
        record = f.read_record('<i4', '<f4')

        n_vis = record[0][0]
        time = record[1][0]

        # Read point data

        Teff = np.empty(n_vis)
        V_proj = np.empty(n_vis)
        A_proj = np.empty(n_vis)
        g = np.empty(n_vis)
        mu = np.empty(n_vis)

        for i in range(n_vis):
            try:
                Teff[i], V_proj[i], A_proj[i], g[i], mu[i] = f.read_reals('<f4')
            except si.FortranEOFError:
                break

        # Initialize the table

        tbl = at.Table({
            'Teff': Teff[:i],
            'V_proj': V_proj[:i],
            'A_proj': A_proj[:i],
            'g': g[:i],
            'mu': mu[:i]},
            units={'K', 'm s^-1', 'm^2', 'm s^-2', ''},
            meta={'n_vis': i, 'time': time}
        )

    # Return the table

    return tbl


def parse_comm_file(file_name, check_against=None):

    r"""Read and parse a BRUCE/KYLIE command file.

    Args:
        file_name (str): Name of file.
        check_against (dict, optional): Grammar rules to check against.

    Returns:
        dict: Nested data structure with commands as the top-level keys. Indexing
        is dict[comm][i][param], where comm is the command, i is the index (to
        handle multiple commands with the same name), and param is the parameter.
    """

    # Read the file into a list of tokens

    tokens = []

    with open(file_name, 'r') as f:

        tokens = []

        comm = None

        # Process each line

        for line in f:

            # Strip leading/trailing whitespace

            line = line.strip()

            # Strip comments

            line = re.sub(r'[!%].*', '', line)

            # Skip if empty

            if not line:
                continue

            # Append to the tokens list

            tokens += line.split()

    # Now parse the token list to build the command/parameter dict

    comms = {}
    
    comm = None
    params = None

    for token in tokens:

        if comm:

            if token[0] == '#':
                raise Exception(f'Misplaced # inside command {comm} when parsing file {file_name}')
            elif token == '{':
                if params:
                    raise Exception(f'Misplaced {{ inside command {comm} when parsing file {file_name}')
                params = {}
            elif token == '}':
                if not params:
                    raise Exception(f'Misplaced }} inside command {comm} when parsing file {file_name}')
                comms[comm] += [params]
                params = None
                comm = None
            else:
                try:
                    param, value = token.split(':')
                    try:
                        value = int(value)
                    except ValueError:
                        try:
                            value = float(value)
                        except ValueError:
                            pass
                    params[param] = value
                except ValueError:
                    params[token] = True

        else:

            if token[0] == '#':
                comm = token[1:]
                if not comm in comms:
                    comms[comm] = []
                params = None
            else:
                raise Exception(f'Missing # when parsing file {file_name}')

    # If necessary, perform checks

    if check_against is not None:

        for comm in check_against:

            n_min, n_max = check_against[comm]

            if n_min is not None:
                if not comm in comms:
                    raise Exception(f'Missing {comm} command when parsing file {file_name}')
                elif len(comms[comm]) < n_min:
                    raise Exception(f'Too few {comm} commands when parsing file {file_name}')

            if n_max is not None:
                if comm in comms:
                    if len(comms[comm]) > n_max:
                        raise Exception(f'Too many {comm} commands when parsing file {file_name}')

    # Return the commands

    return comms


def integrate_flux(tbl, specgrid, lam, d=10, limb_u=None, x_add=None):

    r"""Disk-integrate a flux spectrum.

    Args:
        tbl (astropy.table.Table): Table of BRUCE model data returned
            by :py:func:`read_bruce_model`.
        specgrid (pymsg.SpecGrid): Grid of spectroscopic intensity data.
        lam (numpy.ndarray): Wavelength abscissa (Å).
        d (float, optional): Distance (pc).
        limb_u (float, optional): Linear limb-darkening parameter.
        x_add (dict, optional): Additional photospheric parameters (e.g., [Fe/H])
            passed on to pymsg.

    Returns:
        numpy.ndarray: Spectroscopic flux (erg/cm^2/s/Å) in bins delineated
        by lam; length len(lam)-1.
    """

    # Initialize the flux array

    flux = np.zeros(len(lam)-1)

    # Add in contributions from each surface element

    i = 0

    for row in tbl:

        # Set up photospheric parameters

        if x_add is not None:
            x = x_add | {
                'Teff': row['Teff'],
                'log(g)': np.log10(row['g'])+2
            }
        else:
            x = {
                'Teff': row['Teff'],
                'log(g)': np.log10(row['g'])+2
            }

        # Doppler-shift the wavelength axis

        lam_d = lam/(1-row['V_proj']/ac.c.value)

        # Add flux contribution

        if limb_u is not None:
            intensity = specgrid.flux(x, lam_d)*(1 - limb_u*(1 - row['mu']))/(np.pi*(1 - limb_u/3))
        else:
            intensity = specgrid.intensity(x, row['mu'], lam_d)

        flux += intensity*row['A_proj']/(d*ac.pc.value)**2

        i += 1
        
    # Return the flux

    return flux


def run_kylie(file_name):

    '''Perform a KYLIE run.

    Args:
        file_name (str): Name of KYLIE input file.
    '''
    
    # Read the KYLIE input file

    comms = parse_comm_file(file_name,
                            {'fields': [1,1],
                             'waveband': [None,None],
                             'wavepoint': [None,None],
                             'specgrid': [1,1]})

    # Set up the wavelength abcissae. These is designed to place the
    # central wavelength of each bin at the same locations that KYLIE
    # would place sample points

    lams = []

    if 'waveband' in comms:

        for params in comms['waveband']:

            lam_min = params['start_wavelength']
            lam_max = params['finish_wavelength']
            dlam = params['wavelength_resolution']

            n = int(round((lam_max - lam_min)/dlam)) + 1

            lams += [np.linspace(lam_min-dlam/2, lam_max+dlam/2, n+1)]

    if 'wavepoint' in comms:

        for params in comms['wavepoint']:

            R = 1e6 # Nominal resolving power
            lam_min = params['wavelength']*(1 - 1/(2*R))
            lam_max = params['wavelength']*(1 + 1/(2*R))

            lams += [np.array([lam_min, lam_max])]

    if len(lams) == 0:
        raise Exception('No wavelength points defined')

    lam_min = min([min(lam) for lam in lams])
    lam_max = min([max(lam) for lam in lams])

    # Load the specgrid

    specgrid = pm.SpecGrid(comms['specgrid'][0]['file_name'])

    x_add = {}

    for param in comms['specgrid'][0]:
        if param != 'file_name':
            x_add[param] = comms['specgrid'][0][param]

    # Set caching parameters, based on the wavelength range in lams
    # with a 1000 km/s margin of error

    moe = 1000e3

    specgrid.cache_lam_min = lam_min/(1 + moe/ac.c.value)
    specgrid.cache_lam_max = lam_max/(1 - moe/ac.c.value)

    # If necessary, set the limb-darkening parameter

    if 'limb_u_override' in comms['fields'][0]:
        limb_u = comms['fields'][0]['limb_u_override']
    else:
        limb_u = None

    # Loop over fields

    for i in range(comms['fields'][0]['number_of_fields']):

        # Set up file names

        model_file_name = f'{comms["fields"][0]["dump_filestub"]}{i+1:03d}'
        spec_file_name = f'{model_file_name}.ecsv'

        # Read the model data

        tbl_model = read_bruce_model(model_file_name)

        # Create the spectrum

        spec_lam = np.empty((0))
        spec_flux = np.empty((0))

        for lam in lams:

            flux = integrate_flux(tbl_model, specgrid, lam, limb_u=limb_u, x_add=x_add)

            spec_lam = np.append(spec_lam, 0.5*(lam[1:] + lam[:-1]))
            spec_flux = np.append(spec_flux, integrate_flux(tbl_model, specgrid, lam, limb_u=limb_u, x_add=x_add))

        # Write the spectrum

        tbl_spec = at.Table({'wavelength': spec_lam, 'flux': spec_flux},
                            units=('Angstrom', 'erg cm^-2 s^-1 Angstrom^-1'))

        tbl_spec.write(spec_file_name, overwrite=True)

        print(f'Written file {spec_file_name}')
