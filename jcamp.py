# -*- coding: UTF-8 -*-
# See the LICENSE.rst file for licensing information.

import datetime
from numpy import array, append, arange, logical_not, log10, nan
import re
import pdb

'''
jcamp.py contains functions useful for parsing JCAMP-DX formatted files containing spectral data. The main
function `jcamp_readfile()` formats the input file into a Python dictionary, while `jcamp_calc_xsec()`
converts a given JCAMP-style data dictionary from absorption units to cross-section (m^2).

The bottom of the file contains an example script, so that if the module is run by itself, it will show several
spectra plotted from data in repository folders.
'''

__authors__ = 'Nathan Hagen'
__license__ = 'MIT/X11 License'
__contact__ = 'Nathan Hagen <and.the.light.shattered@gmail.com>'
__all__     = ['jcamp_readfile', 'jcamp_calc_xsec', 'is_float', 'get_value', 'jcamp_parse']
__version__ = '1.2.1'

## In SQZ_digits, '+' or '-' is for PAC, ',' for CSV.
SQZ_digits = {'@':'+0', 'A':'+1', 'B':'+2', 'C':'+3', 'D':'+4', 'E':'+5', 'F':'+6', 'G':'+7', 'H':'+8', 'I':'+9',
              'a':'-1', 'b':'-2', 'c':'-3', 'd':'-4', 'e':'-5', 'f':'-6', 'g':'-7', 'h':'-8', 'i':'-9',
              '+':'+',  '-':'-',  ',':' '}
DIF_digits = {'%': 0, 'J':1,  'K':2,  'L':3,  'M':4,  'N':5,  'O':6,  'P':7,  'Q':8,  'R':9,
              'j':-1, 'k':-2, 'l':-3, 'm':-4, 'n':-5, 'o':-6, 'p':-7, 'q':-8, 'r':-9}
DUP_digits = {'S':1, 'T':2, 'U':3, 'V':4, 'W':5, 'X':6, 'Y':7, 'Z':8, 's':9}

# The specification allows multiple formats for representing LONGDATE.
# See `FRACTIONAL_SECONDS_PATTERN` below for the optional token representing fractional seconds.
# These fractional seconds are removed in advance. Thus `%N` is not referenced in the formats below.
DATE_FORMATS = ["%Y/%m/%d %H:%M:%S %z", "%Y/%m/%d %H:%M:%S", "%Y/%m/%d"]

# The optional token describing the fractional seconds is referenced in the specification as `.SSSS`.
# This number of digits (four) is rather unclear, since the usual presentation of a fraction of
# seconds would contain either 3, 6 or 9 digits.
FRACTIONAL_SECONDS_PATTERN = re.compile(
    r"^\d{4}/\d{2}/\d{2} +\d{2}:\d{2}\d{2}(?P<fractional_seconds>\d{1,9})"
)

##=====================================================================================================
def jcamp_readfile(filename):
    with open(filename, 'rb') as filehandle:
        data = jcamp_read(filehandle)
    data['filename'] = filename
    return(data)

##=====================================================================================================
def _parse_longdate(date_string: str) -> datetime.datetime:
    """parse the "LONGDATE" field according to the JCAMP-DX specification

    raises ValueError in case of problems
    """
    fractional_seconds_match = FRACTIONAL_SECONDS_PATTERN.search(date_string)
    if fractional_seconds_match:
        # Remove the fractional seconds string - it would complicate `strptime`.
        date_string = FRACTIONAL_SECONDS_PATTERN.sub("", date_string)

        # Try to interprete the fractional seconds. The JCAMP specification (v6.00) does not
        # explain, how a string of arbitrary length is supposed to be interpreted.
        # Thus we are just guessing based on the number of digits.
        fraction_seconds_string = fractional_seconds_match.group("fractional_seconds")
        if len(fraction_seconds_string) in {7, 8, 9}:
            # this is probably nanoseconds
            microseconds = int(int(fraction_seconds_string) / 1000)
        elif len(fraction_seconds_string) in {4, 5, 6}:
            microseconds = int(fraction_seconds_string)
        elif len(fraction_seconds_string) in {1, 2, 3}:
            microseconds = 1000 * int(fraction_seconds_string)
        else:
            # We should never end up here.
            raise ValueError("Fractional seconds string could not be parsed: {}".format(fraction_seconds_string))
    else:
        microseconds = 0

    # Parse the date and time.
    for fmt in DATE_FORMATS:
        try:
            parsed = datetime.datetime.strptime(date_string, fmt)
        except ValueError:
            pass
        else:
            # Inject the previously parsed microseconds
            return parsed.replace(microsecond=microseconds)
    else:
        raise ValueError("Failed to parse the date string: {}".format(date_string))

##=====================================================================================================
def jcamp_read(filehandle):
    '''
    Read a JDX-format file, and return a dictionary containing the header info, a 1D numpy vectors `x` for the
    abscissa information (e.g. wavelength or wavenumber) and `y` for the ordinate information (e.g. transmission).

    Parameters
    ----------
    filehandle : str
        The object representing the JCAMP-DX filename to read.

    Returns
    -------
    jcamp_dict : dict
        The dictionary containing the header and data vectors.
    '''

    jcamp_dict = {}
    xstart = []
    xnum = []
    y = []
    x = []
    datastart = False
    is_compound = False
    in_compound_block = False
    compound_block_contents = []
    re_num = re.compile(r'\d+')
    lhs = None
    for line in filehandle:
        ## When parsing compound files, the input is an array of strings, so no need to decode it twice.
        if hasattr(line, 'decode'):
            line = line.decode('utf-8','ignore')

        if not line.strip():
            continue
        if line.startswith('$$'):
            continue

        ## Detect the start of a compound block
        if is_compound and line.upper().startswith('##TITLE'):
            in_compound_block = True
            compound_block_contents = [line]
            continue

        ## If we are reading a compound block, collect lines into an array to be processed by a
        ## recursive call this this function.
        if in_compound_block:
            ## Store this line.
            compound_block_contents.append(line)

            ## Detect the end of the compound block.
            if line.upper().startswith('##END'):
                ## Process the entire block and put it into the children array.
                jcamp_dict['children'].append(jcamp_read(compound_block_contents))
                in_compound_block = False
                compound_block_contents = []
            continue

        ## Lines beginning with '##' are header lines.
        if line.startswith('##'):
            line = line.strip('##')
            (lhs,rhs) = line.split('=', 1)
            lhs = lhs.strip().lower()
            rhs = rhs.strip()
            #continuation = rhs.endswith('=')

            if rhs.isdigit():
                jcamp_dict[lhs] = int(rhs)
            elif is_float(rhs):
                jcamp_dict[lhs] = float(rhs)
            else:
                jcamp_dict[lhs] = rhs

            ## Detect compound files.
            ## See table XI in http://www.jcamp-dx.org/protocols/dxir01.pdf
            if (lhs in {'data type', 'datatype'}) and (rhs.lower() == 'link'):
                is_compound = True
                jcamp_dict['children'] = []

            if (lhs in ('xydata', 'xypoints', 'peak table')):
                ## This is a new data entry, reset x and y.
                x = []
                y = []
                datastart = True
                datatype = rhs
                continue        ## data starts on next line
            elif (lhs == 'end'):
                bounds = [int(i) for i in re_num.findall(rhs)]
                datastart = True
                datatype = bounds
                datalist = []
                continue
            elif lhs == 'longdate':
                try:
                    parsed = _parse_longdate(jcamp_dict[lhs])
                except ValueError:
                    # Keep the original date string.
                    pass
                else:
                    # Replace the string with the datetime object.
                    jcamp_dict[lhs] = parsed
            elif datastart:
                datastart = False
        elif lhs is not None and not datastart:  # multiline entry
            jcamp_dict[lhs] += '\n{}'.format(line.strip())

        if datastart and (datatype == '(X++(Y..Y))'):
            ## If the line does not start with '##' or '$$' then it should be a data line.
            ## The pair of lines below involve regex splitting on floating point numbers and integers. We can't just
            ## split on spaces because JCAMP allows minus signs to replace spaces in the case of negative numbers.
            datavals = jcamp_parse(line)
            xstart.append(float(datavals[0]))
            xnum.append(len(datavals) - 1)
            for dataval in datavals[1:]:
                y.append(float(dataval))
        elif datastart and (('xypoints' in jcamp_dict) or ('xydata' in jcamp_dict)) and (datatype == '(XY..XY)'):
            datavals = [v.strip() for v in re.split(r"[,;\s]", line) if v]  ## be careful not to allow empty strings
            if not all(is_float(datavals)): continue
            datavals = array(datavals)
            x.extend(datavals[0::2])        ## every other data point starting at the zeroth
            y.extend(datavals[1::2])        ## every other data point starting at the first
        elif datastart and ('peak table' in jcamp_dict) and (datatype == '(XY..XY)'):
            datavals = [v.strip() for v in re.split(r"[,;\s]", line) if v]  ## be careful not to allow empty strings
            if not all(is_float(datavals)): continue
            datavals = array(datavals)
            x.extend(datavals[0::2])        ## every other data point starting at the zeroth
            y.extend(datavals[1::2])        ## every other data point starting at the first
        elif datastart and isinstance(datatype,list):
            ## If the line does not start with '##' or '$$' then it should be a data line.
            ## The pair of lines below involve regex splitting on floating point numbers and integers. We can't just
            ## split on spaces because JCAMP allows minus signs to replace spaces in the case of negative numbers.
            datavals = jcamp_parse(line)
            datalist += datavals

    if ('xydata' in jcamp_dict) and (jcamp_dict['xydata'] == '(X++(Y..Y))'):
        ## You got all of the Y-values. Next you need to figure out how to generate the missing X's...
        ## First look for the "lastx" dictionary entry. You will need that one to finish the set.
        xstart.append(jcamp_dict['lastx'])
        x = array([])
        for n in range(len(xnum)-1):
            dx = (xstart[n+1] - xstart[n]) / xnum[n]
            x = append(x, xstart[n]+(dx*arange(xnum[n])))
            #print(n, xstart[n], xstart[n+1], xnum[n], xstart[n]+(dx*arange(xnum[n])))

        ## The last line must be treated separately.
        if (xnum[len(xnum)-1] > 1):
            dx = (jcamp_dict['lastx'] - xstart[len(xnum)-1]) / (xnum[len(xnum)-1] - 1.0)
            x = append(x, xstart[len(xnum)-1]+(dx*arange(xnum[len(xnum)-1])))
            #print(n, xstart[len(xnum)-1]+(dx*arange(xnum[len(xnum)-1])))
        else:
            x = append(x, jcamp_dict['lastx'])

        y = array([float(yval) for yval in y])
    else:
        x = array([float(xval) for xval in x])
        y = array([float(yval) for yval in y])

    ## The "xfactor" and "yfactor" variables contain any scaling information that may need to be applied
    ## to the data. Go ahead and apply them.
    if ('xfactor' in jcamp_dict):
        x = x * jcamp_dict['xfactor']
    if ('yfactor' in jcamp_dict):
        y = y * jcamp_dict['yfactor']
    jcamp_dict['x'] = x
    jcamp_dict['y'] = y

    return(jcamp_dict)

##=====================================================================================================
def jcamp_calc_xsec(jcamp_dict, wavemin=None, wavemax=None, skip_nonquant=True, debug=False):
    '''
    Taking as input a JDX file, extract the spectrum information and transform the absorption spectrum from existing
    units to absorption cross-section.

    This function also corrects for unphysical data (such as negative transmittance values, or transmission above
    1.0), and calculates absorbance if transmittance given. Instead of a return value, the function inserts the
    information into the input dictionary.

    Note that the conversion assumes that the measurements were collected for gas at a temperature of 296K (23 degC).

    Parameters
    ----------
    jcamp_dict : dict
        A JCAMP spectrum dictionary.
    wavemin : float, optional
        The shortest wavelength in the spectrum to limit the calculation to.
    wavemax : float, optional
        The longest wavelength in the spectrum to limit the calculation to.
    skip_nonquant: bool
        If True then return "None" if the spectrum is missing quantitative data. If False, then try to fill in \
        missing quantitative values with defaults.
    '''

    x = jcamp_dict['x']
    y = jcamp_dict['y']

    T = 296.0            ## the temperature (23 degC) used by NIST when collecting spectra
    R = 1.0355E-25       ## the constant for converting data (includes the gas constant)

    ## Note: normally when we convert from wavenumber to wavelength units, the ordinate must be nonuniformly
    ## rescaled in order to compensate. But this is only true if we resample the abscissa to a uniform sampling
    ## grid. In this case here, we keep the sampling grid nonuniform in wavelength space, such that each digital
    ## bin retains its proportionality to energy, which is what we want.
    if (jcamp_dict['xunits'].lower() in ('1/cm', 'cm-1', 'cm^-1')):
        jcamp_dict['wavenumbers'] = array(x)            ## note that array() always performs a copy
        x = 10000.0 / x
        jcamp_dict['wavelengths'] = x
    elif (jcamp_dict['xunits'].lower() in ('micrometers', 'um', 'wavelength (um)')):
        jcamp_dict['wavelengths'] = x
        jcamp_dict['wavenumbers'] = 10000.0 / x
    elif (jcamp_dict['xunits'].lower() in ('nanometers', 'nm', 'wavelength (nm)')):
        x = x / 1000.0
        jcamp_dict['wavelengths'] = x
        jcamp_dict['wavenumbers'] = 10000.0 / x
    else:
        raise ValueError('Don\'t know how to convert the spectrum\'s x units ("' + jcamp_dict['xunits'] + '") to micrometers.')

    ## Correct for any unphysical negative values.
    y[y < 0.0] = 0.0

    ## Make sure "y" refers to absorbance.
    if (jcamp_dict['yunits'].lower() == 'transmittance'):
        ## If in transmittance, then any y > 1.0 are unphysical.
        y[y > 1.0] = 1.0

        ## Convert to absorbance.
        okay = (y > 0.0)
        y[okay] = log10(1.0 / y[okay])
        y[logical_not(okay)] = nan

        jcamp_dict['absorbance'] = y
    elif (jcamp_dict['yunits'].lower() == 'absorbance'):
        pass
    elif (jcamp_dict['yunits'].lower() == '(micromol/mol)-1m-1 (base 10)'):
        jcamp_dict['yunits'] = 'xsec (m^2))'
        jcamp_dict['xsec'] = y / 2.687e19
        return
    else:
        raise ValueError('Don\'t know how to convert the spectrum\'s y units ("' + jcamp_dict['yunits'] + '") to absorbance.')

    ## Determine the effective path length "ell" of the measurement chamber, in meters.
    if ('path length' in jcamp_dict):
        (val,unit) = jcamp_dict['path length'].lower().split()[0:2]
        if (unit == 'cm'):
            ell = float(val) / 100.0
        elif (unit == 'm'):
            ell = float(val)
        elif (unit == 'mm'):
            ell = float(val) / 1000.0
        else:
            ell = 0.1
    else:
        if skip_nonquant: return({'info':None, 'x':None, 'xsec':None, 'y':None})
        ell = 0.1
        if debug: print('Path length variable not found. Using 0.1m as a default ...')

    assert(len(x) == len(y))

    if ('npoints' in jcamp_dict):
        if (len(x) != jcamp_dict['npoints']):
            npts_retrieved = str(len(x))
            msg = '"' + jcamp_dict['title'] + '": Number of data points retrieved (' + npts_retrieved + \
                  ') does not equal the expected length (npoints = ' + str(jcamp_dict['npoints']) + ')!'
            raise ValueError(msg)

    ## For each gas, manually define the pressure "p" at which the measurement was taken (in units of mmHg).
    ## These values are obtained from the NIST Infrared spectrum database, which for some reason did not
    ## put the partial pressure information into the header.
    if ('partial_pressure' in jcamp_dict):
        p = float(jcamp_dict['partial_pressure'].split()[0])
        p_units = jcamp_dict['partial_pressure'].split()[1]
        if (p_units.lower() == 'mmhg'):
            pass
        elif (p_units.lower() == 'ppm'):
            p = p * 759.8 * 1.0E-6       ## scale PPM units at atmospheric pressure to partial pressure in mmHg
    else:
        if debug: print('Partial pressure variable value for ' + jcamp_dict['title'] + ' is missing. Using the default p = 150.0 mmHg ...')
        if skip_nonquant: return({'info':None, 'x':None, 'xsec':None, 'y':None})
        p = 150.0

    ## Convert the absorbance units to cross-section in meters squared per molecule.
    xsec = y * T * R / (p * ell)

    ## Add the "xsec" values to the data dictionary.
    jcamp_dict['xsec'] = xsec

    return

##=====================================================================================================
def is_float(s):
    '''
    Test if a string, or list of strings, contains a numeric value(s).

    Parameters
    ----------
    s : str, or list of str
        The string or list of strings to test.

    Returns
    -------
    is_float_bool : bool or list of bool
        A single boolean or list of boolean values indicating whether each input can be converted into a float.
    '''

    if isinstance(s,tuple) or isinstance(s,list):
        if not all(isinstance(i, str) for i in s):
            raise TypeError("Input {} is not a list of strings".format(s))
        if (len(s) == 0):
            raise ValueError('Input {} is empty'.format(s))
        else:
            bool = list(True for i in range(0,len(s)))
            for i in range(0,len(s)):
                try:
                    float(s[i])
                except ValueError:
                    bool[i] = False
        return(bool)
    else:
        if not isinstance(s, str):
            raise TypeError("Input '%s' is not a string" % (s))

        try:
            float(s)
            return(True)
        except ValueError:
            return(False)

##=====================================================================================================
def get_value(num, is_dif, vals):
    n = float(num)
    if is_dif:
        lastval = vals[-1]
        val = n + lastval
    else:
        val = n

    return(val)

##=====================================================================================================
def jcamp_parse(line):
    line = line.strip()

    datavals = []
    num = ""

    ## Convert whitespace into single space by splitting the string then re-assembling with single spaces.
    line = ' '.join(line.split())

    ## If there are any coded digits, then replace the codes with the appropriate numbers.
    ## 'DUP_digits' are characters that represent how many times the previous character should be replicated.
    ## 'DIF_digits' represent ...?
    ## 'SQZ_digits' represent ...?
    DUP_set = set(DUP_digits)

    if any(c in DUP_set for c in line):
        ## Split the line into individual characters so that you can check for coded characters one-by-one.
        newline = ''
        for (i,c) in enumerate(line):
            if (c in DUP_digits):
                prev_c = line[i-1]
                mul = DUP_digits[c]
                newline += prev_c * (mul-1)
            else:
                mul = ''
                newline += c
        line = "".join(newline)

    DIF = False
    for c in line:
        if c.isdigit() or (c == "."):
            num += c
        elif (c == ' '):
            DIF = False
            if num:
                n = get_value(num, DIF, datavals)
                datavals.append(n)
            num = ''
        elif (c in SQZ_digits):
            DIF = False
            if num:
                n = get_value(num, DIF, datavals)
                datavals.append(n)
            num = SQZ_digits[c]
        elif (c in DIF_digits):
            if num:
                n = get_value(num, DIF, datavals)
                datavals.append(n)
            DIF = True
            num = str(DIF_digits[c])
        else:
            raise Exception("Unknown character (%s) encountered while parsing data" % c)

    if num:
        n = get_value(num, DIF, datavals)
        datavals.append(n)

    return(datavals)

## =================================================================================================
## =================================================================================================

if (__name__ == '__main__'):
    import matplotlib.pyplot as plt
    filename = './data/infrared_spectra/ethylene.jdx'
    jcamp_dict = jcamp_readfile(filename)
    plt.plot(jcamp_dict['x'], jcamp_dict['y'])
    plt.title(filename)
    plt.xlabel(jcamp_dict['xunits'])
    plt.ylabel(jcamp_dict['yunits'])

    jcamp_calc_xsec(jcamp_dict, skip_nonquant=False, debug=False)
    plt.figure()
    plt.plot(jcamp_dict['wavelengths'], jcamp_dict['xsec'])
    plt.title(filename)
    plt.xlabel('wavelength (um)')
    plt.ylabel('absorption cross-section (m^2)')

    filename = './data/uvvis_spectra/toluene.jdx'
    plt.figure()
    jcamp_dict = jcamp_readfile(filename)
    plt.plot(jcamp_dict['x'], jcamp_dict['y'], 'r-')
    plt.title(filename)
    plt.xlabel(jcamp_dict['xunits'])
    plt.ylabel(jcamp_dict['yunits'])

    filename = './data/mass_spectra/ethanol_ms.jdx'
    jcamp_dict = jcamp_readfile(filename)
    plt.figure()
    for n in arange(len(jcamp_dict['x'])):
        plt.plot((jcamp_dict['x'][n],jcamp_dict['x'][n]), (0.0, jcamp_dict['y'][n]), 'm-', linewidth=2.0)
    plt.title(filename)
    plt.xlabel(jcamp_dict['xunits'])
    plt.ylabel(jcamp_dict['yunits'])

    filename = './data/raman_spectra/tannic_acid.jdx'
    jcamp_dict = jcamp_readfile(filename)
    plt.figure()
    plt.plot(jcamp_dict['x'], jcamp_dict['y'], 'k-')
    plt.title(filename)
    plt.xlabel(jcamp_dict['xunits'])
    plt.ylabel(jcamp_dict['yunits'])

    filename = './data/neutron_scattering_spectra/emodine.jdx'
    jcamp_dict = jcamp_readfile(filename)
    plt.figure()
    plt.plot(jcamp_dict['x'], jcamp_dict['y'], 'k-')
    plt.title(filename)
    plt.xlabel(jcamp_dict['xunits'])
    plt.ylabel(jcamp_dict['yunits'])

    filename = './data/infrared_spectra/example_compound_file.jdx'
    jcamp_dict = jcamp_readfile(filename)
    plt.figure()
    for c in jcamp_dict['children']:
        plt.plot(c['x'], c['y'])
    plt.xlabel(jcamp_dict['children'][0]['xunits'])  ## assume all blocks have the same units
    plt.ylabel(jcamp_dict['children'][0]['yunits'])
    plt.title(filename)

    filename = './data/infrared_spectra/example_multiline_datasets.jdx'
    jcamp_dict = jcamp_readfile(filename)
    plt.figure()
    plt.plot(jcamp_dict['x'], jcamp_dict['y'])
    plt.title(filename)
    plt.xlabel(jcamp_dict['xunits'])
    plt.ylabel(jcamp_dict['yunits'])
    print(jcamp_dict['comments'])

    plt.show()
