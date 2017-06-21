import unittest

import numpy
import os
from jcamp import JCAMP_calc_xsec, JCAMP_reader, jcamp_parse


class TestJcamp(unittest.TestCase):

    def assertAlmostEqualRelative(self, val1, val2, tolerance=1e-2):
        diff = float(abs(val1 - val2))
        if val1 and diff / val1 > tolerance and diff / val2 > tolerance:
            raise Exception("%s is not close enough from %s" % (val1, val2))
        return True

    def test_read_IR(self):
        for root, dirs, files in os.walk("./data/infrared_spectra"):
            for filename in files:
                full_filename = os.path.join(root, filename)
                jcamp_dict = JCAMP_reader(full_filename)
                self.assertIsInstance(jcamp_dict['x'], numpy.ndarray)
                self.assertIsInstance(jcamp_dict['y'], numpy.ndarray)
                self.assertEqual(len(jcamp_dict['x']), len(jcamp_dict['y']))
                self.assertEqual(len(jcamp_dict['x']), jcamp_dict['npoints'])
                self.assertAlmostEqualRelative(min(jcamp_dict['x']), jcamp_dict['minx'])
                self.assertAlmostEqualRelative(min(jcamp_dict['y']), jcamp_dict['miny'])
                self.assertAlmostEqualRelative(max(jcamp_dict['x']), jcamp_dict['maxx'])
                self.assertAlmostEqualRelative(max(jcamp_dict['y']), jcamp_dict['maxy'])

    def test_read_raman(self):
        for root, dirs, files in os.walk("./data/raman_spectra"):
            for filename in files:
                full_filename = os.path.join(root, filename)
                jcamp_dict = JCAMP_reader(full_filename)
                self.assertIsInstance(jcamp_dict['x'], numpy.ndarray)
                self.assertIsInstance(jcamp_dict['y'], numpy.ndarray)
                self.assertEqual(len(jcamp_dict['x']), len(jcamp_dict['y']))
                self.assertEqual(len(jcamp_dict['x']), jcamp_dict['npoints'])
                self.assertAlmostEqual(min(jcamp_dict['x']), jcamp_dict['minx'])
                self.assertAlmostEqual(min(jcamp_dict['y']), jcamp_dict['miny'])
                self.assertAlmostEqual(max(jcamp_dict['x']), jcamp_dict['maxx'])
                self.assertAlmostEqual(max(jcamp_dict['y']), jcamp_dict['maxy'])

    def test_read_uv(self):
        for root, dirs, files in os.walk("./data/uvvus_spectra"):
            for filename in files:
                full_filename = os.path.join(root, filename)
                jcamp_dict = JCAMP_reader(full_filename)
                self.assertIsInstance(jcamp_dict['x'], numpy.ndarray)
                self.assertIsInstance(jcamp_dict['y'], numpy.ndarray)
                self.assertEqual(len(jcamp_dict['x']), len(jcamp_dict['y']))
                self.assertEqual(len(jcamp_dict['x']), jcamp_dict['npoints'])
                self.assertAlmostEqual(min(jcamp_dict['x']), jcamp_dict['minx'])
                self.assertAlmostEqual(min(jcamp_dict['y']), jcamp_dict['miny'])
                self.assertAlmostEqual(max(jcamp_dict['x']), jcamp_dict['maxx'])
                self.assertAlmostEqual(max(jcamp_dict['y']), jcamp_dict['maxy'])

    def test_line_parse(self):
        # tests from http://wwwchem.uwimona.edu.jm/software/jcampdx.html
        line = "99 98 97 96 98 93"
        self.assertEquals(jcamp_parse(line), [99, 98, 97, 96, 98, 93])

        line = "99,98,97,96,98,93"
        self.assertEquals(jcamp_parse(line), [99, 98, 97, 96, 98, 93])

        line = "99+98+97+96+98+93"
        self.assertEquals(jcamp_parse(line), [99, 98, 97, 96, 98, 93])

        line = "99I8I7I6I8I3"
        self.assertEquals(jcamp_parse(line), [99, 98, 97, 96, 98, 93])

        line = "99jjjKn"
        self.assertEquals(jcamp_parse(line), [99, 98, 97, 96, 98, 93])

        line = "99jUKn"
        self.assertEquals(jcamp_parse(line), [99, 98, 97, 96, 98, 93])

    def test_spec_line_parse(self):
        # Tests from http://www.jcamp-dx.org/drafts/JCAMP6_2b%20Draft.pdf
        line = "1000 2000 2001 2002 2003 2003 2003"
        self.assertEquals(jcamp_parse(line), [1000, 2000, 2001, 2002, 2003, 2003, 2003])

        line = "+1000+2000+2001+2002+2003+2003+2003"
        self.assertEquals(jcamp_parse(line), [1000, 2000, 2001, 2002, 2003, 2003, 2003])

        line = "A000B000B001B002B003B003B003"
        self.assertEquals(jcamp_parse(line), [1000, 2000, 2001, 2002, 2003, 2003, 2003])

        line = "A000J000JJJ%%"
        self.assertEquals(jcamp_parse(line), [1000, 2000, 2001, 2002, 2003, 2003, 2003])

        # Probably a bug in the example in the spec.
        line = "A000J000JU%%"
        self.assertEquals(jcamp_parse(line), [1000, 2000, 2001, 2002, 2003, 2003, 2003])

    def test_jcamp_reader_dict(self):
        filename = './data/infrared_spectra/methane.jdx'
        jcamp_dict = JCAMP_reader(filename)
        self.assertEqual(jcamp_dict['origin'], 'DOW CHEMICAL COMPANY')
        self.assertEqual(jcamp_dict['deltax'], 0.935748)
        self.assertEqual(jcamp_dict['sampling procedure'], 'TRANSMISSION')
        self.assertEqual(jcamp_dict['firstx'], 449.47)
        self.assertEqual(jcamp_dict['firsty'], 0.953)
        self.assertEqual(jcamp_dict['data type'], 'INFRARED SPECTRUM')
        self.assertEqual(jcamp_dict['owner'], 'COBLENTZ SOCIETY')
        self.assertEqual(jcamp_dict['maxx'], 3801.32)
        self.assertEqual(jcamp_dict['maxy'], 1.03697)
        self.assertEqual(jcamp_dict['end'], [])
        self.assertEqual(jcamp_dict['title'], 'METHANE')
        self.assertEqual(jcamp_dict['lastx'], 3801.32)
        self.assertEqual(jcamp_dict['state'], 'GAS')
        self.assertEqual(jcamp_dict['yunits'], 'TRANSMITTANCE')
        self.assertEqual(jcamp_dict['cas registry no'], '74-82-8')
        self.assertEqual(jcamp_dict['spectrometer/data system'],
                         'DOW KBr FOREPRISM')
        self.assertEqual(jcamp_dict['$nist image'], 'cob8873')
        self.assertEqual(jcamp_dict['path length'], '5 CM')
        self.assertEqual(jcamp_dict['data processing'],
                         'DIGITIZED BY NIST FROM HARD COPY (FROM TWO SEGMENTS)'
                         )
        self.assertEqual(jcamp_dict['date'], 1964)
        self.assertEqual(jcamp_dict['molform'], 'C H4')
        self.assertEqual(jcamp_dict['instrument parameters'],
                         'GRATING CHANGED AT 5.0, 7.5, 15.0 MICRON')
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
        jcamp_dict = JCAMP_reader(filename)
        JCAMP_calc_xsec(jcamp_dict)
        self.assertTrue('xsec' in jcamp_dict)


if __name__ == '__main__':
    unittest.main()
