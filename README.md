jcamp
=====

A set of Python utilities for reading JCAMP-DX files.

The ``JCAMP_reader()`` function takes a filename as input, and returns a dictionary containing the data found in the file. Specifically, the keys contained in the dictionary are: (1) the field names found in the file's header, with values being int- or float-type if the corresponding field is a numerical type, or a string-type otherwise. (2) two arrays ``x`` and ``y``, giving the scaled values of the data points (scaled according to the ``xfactor`` and ``yfactor`` fields in the header, if they exist. The units of ``x`` and ``y`` are whatever are indicated in the header fields ``xunits`` and ``yunits``, if these exist.

The ``JCAMP_calc_xsec()`` function is intended to takes as input the result of the ``JCAMP_reader()`` function and to convert the ``x`` data to wavelength in microns, and the ``y`` data to cross-section in units of m^2 for gas concentration of 1ppm at standard atmospheric pressure and temperature, across a path length of 1 meter. The ``JCAMP_calc_xsec()`` function takes as input the data dictionary ``jcamp_dict``, and manipulates that dictionary directly without having a separate return value. Changes to the dictionary may including adding the fields::

    wavelengths: the array of wavelength values (in microns) for each data point
    wavenumbers: the array of wavenumber values (in cm^-1) for each data point
    xsec: the array of cross-section values (in units of m^2 at 1ppm.m)

and modifying the fields::

    xunits: micron
    yunits: m^2 at 1ppm.m

The optional arguments ``wavemin``, ``wavemax`` are used if the user wishes to truncate the data to only a desired spectral range. For example, setting ``wavemin=8.0`` and ``wavemax=12.0`` means that the returned data arrays will only contain data corresponding to those wavlengths. If the ``skip_nonquant`` optional input argument is used, then any input spectrum that does not have the complete ``path_length`` and ``partial_pressure`` fields in its dictionary will be passed without modification. (That is, no conversion to quantitative cross-section will be attempted.) If this option is set to True, then if the code finds missing data, it will attempt to generate a quantitative cross-section by guessing the missing values. Based upon NIST's infrared database, typical values for guessing here are ``partial_pressure = 150.0 mmHg`` and ``path length = 0.1 m``.

You can view a notebook demo in the doc folder to see how you can produce a series of plots showing spectra.

jcamp files
-----------

The repository comes with four folders containing JCAMP-format files: ``infrared_spectra/``, ``mass_spectra/``, ``raman_spectra/``, and ``uvvis_spectra``. These were downloaded from freely-available internet databases, and can be used as example format files.
