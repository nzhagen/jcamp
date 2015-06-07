import unittest
from jcamp import JCAMP_reader
import numpy


class TestJcampReader(unittest.TestCase):
    def setUp(self):
        filename = './infrared_spectra/methane.jdx'
        self.jcamp_dict = JCAMP_reader(filename)

    def test_jcamp_reader_dict(self):
        jcamp_dict = self.jcamp_dict
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

    def test_jcamp_reader_data(self):
        jcamp_dict = self.jcamp_dict
        self.assertIsInstance(jcamp_dict['x'], numpy.ndarray)
        self.assertIsInstance(jcamp_dict['y'], numpy.ndarray)
        self.assertEqual(len(jcamp_dict['x']), len(jcamp_dict['y']))
        self.assertEqual(len(jcamp_dict['x']), jcamp_dict['npoints'])
        self.assertAlmostEqual(min(jcamp_dict['x']), jcamp_dict['minx'])
        self.assertAlmostEqual(min(jcamp_dict['y']), jcamp_dict['miny'])
        self.assertAlmostEqual(max(jcamp_dict['x']), jcamp_dict['maxx'])
        self.assertAlmostEqual(max(jcamp_dict['y']), jcamp_dict['maxy'], places=4)


if __name__ == '__main__':
    unittest.main()
