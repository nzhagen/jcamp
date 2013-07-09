from numpy import array, linspace, amin, amax, alen

##=====================================================================================================
def JCAMP_reader(filename):
    '''
    Read a JDX-format file, and return a dictionary containing the header info, a 1D numpy vectors `x` for
    the abscissa information (e.g. wavelength or wavenumber) and `y` for the ordinate information (e.g.
    transmission).

    Parameters
    ----------
    filename : str
        The JCAMP-DX filename to read.

    Returns
    -------
    jcamp_dict : dict
        The dictionary containing the header and data vectors.
    '''

    f = open(filename, 'r')
    jcamp_dict = {}
    xstart = []
    xnum = []
    y = []

    for line in f:
        if not line.strip(): continue
        if line.startswith('$$'): continue

        ## Lines beginning with '#' are header lines.
        if line.startswith('#'):
            line = line.strip('#')
            #if line.endswith('='):
            #    lhs = line
            #    rhs = ''
            #else:
            (lhs,rhs) = line.split('=', 1)
            lhs = lhs.strip()
            rhs = rhs.strip()

            if rhs.isdigit():
                jcamp_dict[lhs] = int(rhs)
            elif is_float(rhs):
                jcamp_dict[lhs] = float(rhs)
            else:
                jcamp_dict[lhs] = rhs
        else:
            ## If the line does not start with '##' or '$$' then it should be a data line.
            ## To make sure, check that all split elements are floats.
            datavals = line.split(' ')
            if all(is_float(datavals)):
                xstart.append(float(datavals[0]))
                xnum.append(len(datavals)-1)
                for dataval in datavals[1:]:
                    y.append(float(dataval))

    ## You got all of the Y-values. Next you need to figure out how to generate the missing X's...
    ## First look for the "lastx" dictionary entry. You will need that one to finish the set.
    xstart.append(jcamp_dict['lastx'])
    x = array([])
    y = array(y)
    for n in arange(len(xnum)):
        x = append(x, linspace(xstart[n],xstart[n+1],xnum[n]))

    jcamp_dict['x'] = x
    jcamp_dict['y'] = y

    return(jcamp_dict)

##=====================================================================================================
def JCAMP_calc_xsec(jcamp_dict, wavemin=None, wavemax=None, skip_nonquant=True, debug=False):
    '''
    Taking as input a JDX file, extract the spectrum information and transform the absorption spectrum
    from existing units to absorption cross-section.

    This function also corrects for unphysical data (such as negative transmission values, or
    transmission above 1.0).

    Parameters
    ----------
    jcamp_dict : dict
        A JCAMP spectrum dictionary.
    wavemin : float, optional
        ?
    wavemax : float, optional
        ?
    skip_nonquant: bool
        If True then return "None" if the spectrum is missing quantitative data. If False, then try \
        to fill in missing quantitative values with defaults.

    Returns
    -------
    jcamp_dict : dict
        The input `jcamp_dict` data with the cross-section data `xsec`, `wavelengths`, and `wavenumbers` \
        added as new entries.
    '''

    x = jcamp_dict['x']
    y = jcamp_dict['y']

    ppm = 1.0E-6        ## for conversion to ppm units
    amagat = 2.687E+25  ## number of molecules per m^3 at std temp and pressure
    T = 288.0           ## standard temperature in Kelvin
    R = 2.78            ## the gas constant in (m^3 * mmHg) / (amg * K)

    ## Note: normally when we convert from wavenumber to wavelength units, the ordinate must be nonuniformly
    ## rescaled in order to compensate. But this is only true if we resample the abscissa to a uniform sampling
    ## grid. In this case here, we keep the sampling grid nonuniform in wavelength space, such that each digital
    ## bin retains its proportionality to energy, which is what we want.
    if (jcamp_dict['xunits'] == '1/cm') or (jcamp_dict['xunits'] == 'cm-1'):
        jcamp_dict['wavenumbers'] = array(x)            ## note that array() always performs a copy
        x = 10000.0 / x
        jcamp_dict['wavelengths'] = x
    elif (jcamp_dict['xunits'] == 'micrometers'):
        jcamp_dict['wavelengths'] = x
        jcamp_dict['wavenumbers'] = 10000.0 / x
    elif (jcamp_dict['xunits'] == 'nanometers'):
        x = x * 1000.0
        jcamp_dict['wavelengths'] = x
        jcamp_dict['wavenumbers'] = 10000.0 / x
    else:
        raise ValueError('Don\'t know how to convert the spectrum\'s x units ("' + jcamp_dict['xunits'] + '") to micrometers.')

    ## Make sure "y" refers to absorbance.
    if (jcamp_dict['yunits'] == 'transmittance'):
        y = 1.0 - y
    elif (jcamp_dict['yunits'] == 'absorbance'):
        pass
    else:
        raise ValueError('Don\'t know how to convert the spectrum\'s y units ("' + jcamp_dict['yunits'] + '") to absorbance.')

    ## Correct for any unphysical data.
    if any(y < 0.0):
        y[y < 0.0] = 0.0
    if any(y > 1.0):
        y[y > 1.0] = 1.0

    ## Determine the effective path length "ell" of the measurement chamber, in meters.
    if ('path length' in jcamp_dict):
        (val,unit) = jcamp_dict['path length'].split()[0:2]
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

    assert(alen(x) == alen(y))

    if ('xfactor' in jcamp_dict): x = x * jcamp_dict['xfactor']
    if ('yfactor' in jcamp_dict): y = y * jcamp_dict['yfactor']
    if ('npoints' in jcamp_dict):
        if (alen(x) != jcamp_dict['npoints']):
            npts_retrieved = str(alen(x))
            msg = '"' + jcamp_dict['title'] + '": Number of data points retrieved (' + npts_retrieved + \
                  ') does not equal the expected length (npoints = ' + str(jcamp_dict['npoints']) + ')!'
            raise ValueError(msg)

    ## For each gas, manually define the pressure "p" at which the measurement was taken (in units of mmHg).
    ## These values are obtained from the NIST Infrared spectrum database, which for some reason did not
    ## put the partial pressure information into the header.
    partial_pressures = {'carbon dioxide':200.0,
                         'carbon monoxide':400.0,
                         'hydrogen sulfide':600.0,
                         'methane':150.0,
                         'ethane':150.0,
                         'propane':150.0,
                         'n-pentane':6.0,
                         'propane, 2-methyl-':200.0,        ## also called "iso-butane"
                         'iso-butane':200,
                         'butane, 2-methyl-':100.0,         ## also called "iso-pentane"
                         'iso-pentane':100.0,
                         'ammonia':50.0,
                         'benzene':70.0,
                         'butadiene':100.0,
                         'chlorobenzene':10.0,
                         '1,2-dichloroethane':50.0,
                         'ethanol':30.0,
                         'methanol':70.0,
                         'propylene':150.0,
                         'propylene oxide':100.0,
                         'toluene':20.0,
                         'vinyl chloride':6.0,
                         'p-xylene':5.0,
                         'm-xylene':30.0,
                         'sulfur dioxide':100.0,
                         'butane':100.0,
                         'sulfur hexafluoride':600.0,
                         'ethylene':759.8,
                         'ozone':759.8*40.0/1.0E6}

    if (jcamp_dict['title'] in partial_pressures):
        p = partial_pressures[jcamp_dict['title']]
    else:
        if debug: print('No pressure "p" value entry for ' + jcamp_dict['title'] + '. Using the default p = 150.0 ...')
        if skip_nonquant: return([None]*4)
        p = 150.0

    ## Convert the absorbance units to cross-section in meters squared, for a gas at 1ppm at std atmospheric
    ## temperature and pressure.
    xsec = y * ppm * R * T / (p * ell)

    ## Update the "x" and "y" values and add the "xsec" values.
    jcamp_dict['xsec'] = xsec

    return(jcamp_dict)

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
        if not all(isinstance(i,str) for i in s): raise TypeError("Input 's' is not a list of strings")
        if len(s) == 0:
            try:
                temp = float(i)
            except ValueError:
                return(False)
        else:
            bool = list(True for i in xrange(0,len(s)))
            for i in xrange(0,len(s)):
                try:
                    temp = float(s[i])
                except ValueError:
                    bool[i] = False
        return(bool)
    else:
        if not isinstance(s,str): raise TypeError("Input 's' is not a string")
        try:
            temp = float(s)
            return(True)
        except ValueError:
            return(False)
