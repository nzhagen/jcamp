import unittest
from numpy import amin, amax, ndarray
import os
from jcamp import jcamp_calc_xsec, jcamp_readfile, jcamp_parse
import pdb

class TestJcamp(unittest.TestCase):
    def assertAlmostEqualRelative(self, val1, val2, tolerance=1.0e-2):
        diff = float(abs(val1 - val2))
        if val1 and (diff / val1 > tolerance) and (diff / val2 > tolerance):
            raise Exception("%s is not close enough to %s" % (val1, val2))
        return(True)

    def test_xy_minmax(self, testdict):
        ## Note that 'mass' files seem to be more complex, and current parsing fails on some assumptions.
        self.assertIsInstance(testdict['x'], ndarray)
        self.assertIsInstance(testdict['y'], ndarray)
        self.assertEqual(len(testdict['x']), len(testdict['y']))
        if (len(testdict['x']) > 0) and ('minx' in testdict):
            self.assertEqual(len(testdict['x']), testdict['npoints'])
            self.assertAlmostEqualRelative(amin(testdict['x']), testdict['minx'])
            self.assertAlmostEqualRelative(amin(testdict['y']), testdict['miny'])
            self.assertAlmostEqualRelative(amax(testdict['x']), testdict['maxx'])
            self.assertAlmostEqualRelative(amax(testdict['y']), testdict['maxy'])
        elif ('children' in testdict):
            for child in testdict['children']:
                if ('minx' in testdict):
                    self.assertAlmostEqualRelative(amin(child['x']), child['minx'])
                    self.assertAlmostEqualRelative(amin(child['y']), child['miny'])
                    self.assertAlmostEqualRelative(amax(child['x']), child['maxx'])
                    self.assertAlmostEqualRelative(amax(child['y']), child['maxy'])

    def test_read_IR(self):
        for root, dirs, files in os.walk("./data/infrared_spectra"):
            for filename in files:
                full_filename = os.path.join(root, filename)
                jcamp_dict = jcamp_readfile(full_filename)
                self.test_xy_minmax(jcamp_dict)

    def test_read_raman(self):
        for root, dirs, files in os.walk("./data/raman_spectra"):
            for filename in files:
                full_filename = os.path.join(root, filename)
                jcamp_dict = jcamp_readfile(full_filename)
                self.test_xy_minmax(jcamp_dict)

    def test_read_hnmr_spectra(self):
        for root, dirs, files in os.walk("./data/hnmr_spectra"):
            for filename in files:
                full_filename = os.path.join(root, filename)
                jcamp_dict = jcamp_readfile(full_filename)
                self.test_xy_minmax(jcamp_dict)

    def test_read_uv(self):
        for root, dirs, files in os.walk("./data/uvvus_spectra"):
            for filename in files:
                full_filename = os.path.join(root, filename)
                jcamp_dict = jcamp_readfile(full_filename)
                self.test_xy_minmax(jcamp_dict)

    def test_read_mass(self):
        for root, dirs, files in os.walk("./data/mass_spectra"):
            for filename in files:
                full_filename = os.path.join(root, filename)
                jcamp_dict = jcamp_readfile(full_filename)
                self.test_xy_minmax(jcamp_dict)

    def test_line_parse(self):
        ## Tests from http://wwwchem.uwimona.edu.jm/software/jcampdx.html
        line = "99 98 97 96 98 93"
        self.assertEqual(jcamp_parse(line), [99, 98, 97, 96, 98, 93])
        line = "99,98,97,96,98,93"
        self.assertEqual(jcamp_parse(line), [99, 98, 97, 96, 98, 93])
        line = "99+98+97+96+98+93"
        self.assertEqual(jcamp_parse(line), [99, 98, 97, 96, 98, 93])
        line = "99I8I7I6I8I3"
        self.assertEqual(jcamp_parse(line), [99, 98, 97, 96, 98, 93])
        line = "99jjjKn"
        self.assertEqual(jcamp_parse(line), [99, 98, 97, 96, 98, 93])
        line = "99jUKn"
        self.assertEqual(jcamp_parse(line), [99, 98, 97, 96, 98, 93])

    def test_spec_line_parse(self):
        # Tests from http://www.jcamp-dx.org/drafts/JCAMP6_2b%20Draft.pdf
        line = "1000 2000 2001 2002 2003 2003 2003"
        self.assertEqual(jcamp_parse(line), [1000, 2000, 2001, 2002, 2003, 2003, 2003])
        line = "+1000+2000+2001+2002+2003+2003+2003"
        self.assertEqual(jcamp_parse(line), [1000, 2000, 2001, 2002, 2003, 2003, 2003])
        line = "A000B000B001B002B003B003B003"
        self.assertEqual(jcamp_parse(line), [1000, 2000, 2001, 2002, 2003, 2003, 2003])
        line = "A000J000JJJ%%"
        self.assertEqual(jcamp_parse(line), [1000, 2000, 2001, 2002, 2003, 2003, 2003])
        line = "A000J000JU%%"
        self.assertEqual(jcamp_parse(line), [1000, 2000, 2001, 2002, 2003, 2003, 2003])

    def test_jcamp_readfile_dict(self):
        filename = './data/infrared_spectra/methane.jdx'
        jcamp_dict = jcamp_readfile(filename)
        self.assertEqual(jcamp_dict['origin'], 'DOW CHEMICAL COMPANY')
        self.assertEqual(jcamp_dict['deltax'], 0.935748)
        self.assertEqual(jcamp_dict['sampling procedure'], 'TRANSMISSION')
        self.assertEqual(jcamp_dict['firstx'], 449.47)
        self.assertEqual(jcamp_dict['firsty'], 0.953)
        self.assertEqual(jcamp_dict['data type'], 'INFRARED SPECTRUM')
        self.assertEqual(jcamp_dict['owner'], 'COBLENTZ SOCIETY')
        self.assertEqual(jcamp_dict['maxx'], 3801.32)
        self.assertEqual(jcamp_dict['maxy'], 1.03697)
        self.assertEqual(jcamp_dict['end'], '')
        self.assertEqual(jcamp_dict['title'], 'METHANE')
        self.assertEqual(jcamp_dict['lastx'], 3801.32)
        self.assertEqual(jcamp_dict['state'], 'GAS')
        self.assertEqual(jcamp_dict['yunits'], 'TRANSMITTANCE')
        self.assertEqual(jcamp_dict['cas registry no'], '74-82-8')
        self.assertEqual(jcamp_dict['spectrometer/data system'], 'DOW KBr FOREPRISM')
        self.assertEqual(jcamp_dict['$nist image'], 'cob8873')
        self.assertEqual(jcamp_dict['path length'], '5 CM')
        self.assertEqual(jcamp_dict['data processing'], 'DIGITIZED BY NIST FROM HARD COPY (FROM TWO SEGMENTS)')
        self.assertEqual(jcamp_dict['date'], 1964)
        self.assertEqual(jcamp_dict['molform'], 'C H4')
        self.assertEqual(jcamp_dict['instrument parameters'], 'GRATING CHANGED AT 5.0, 7.5, 15.0 MICRON')
        self.assertEqual(jcamp_dict['class'], 'COBLENTZ')
        self.assertEqual(jcamp_dict['xfactor'], 1.0)
        self.assertEqual(jcamp_dict['partial_pressure'], '150 mmHg')
        self.assertEqual(jcamp_dict['minx'], 449.47)
        self.assertEqual(jcamp_dict['xunits'], '1/CM')
        self.assertEqual(jcamp_dict['source reference'], 'COBLENTZ NO. 8873')
        self.assertEqual(jcamp_dict['miny'], 0.028)
        self.assertEqual(jcamp_dict['jcamp-dx'], 4.24)
        self.assertEqual(jcamp_dict['xydata'], '(X++(Y..Y))')
        self.assertEqual(jcamp_dict['yfactor'], 1)
        self.assertEqual(jcamp_dict['resolution'], 4)
        self.assertEqual(jcamp_dict['$nist source'], 'COBLENTZ')

    def test_jcamp_calc_xsec(self):
        filename = './data/infrared_spectra/methane.jdx'
        jcamp_dict = jcamp_readfile(filename)
        jcamp_calc_xsec(jcamp_dict)
        self.assertTrue('xsec' in jcamp_dict)


if __name__ == '__main__':
    #testobj = TestJcamp()
    unittest.main()
    #testobj.test_read_mass()

