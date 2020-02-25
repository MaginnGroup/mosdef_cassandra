import pytest
import numpy as np

from mosdef_cassandra.tests.base_test import BaseTest
from mosdef_cassandra.utils.convert_box import convert_to_boxmatrix

class TestConvertBox(BaseTest):
    def test_invalid_shape(self):
        box = [10., 5., 5.]
        with pytest.raises(ValueError, match=r"Input must be provided as"):
            box_matrix = convert_to_boxmatrix(box)

    def test_invalid_lengths(self):
        box = [10., 5., -5., 90., 90., 90.]
        with pytest.raises(ValueError, match=r"All box lengths and angles"):
            box_matrix = convert_to_boxmatrix(box)

    def test_invalid_angles(self):
        box = [5., 5., 5., -20., 90., 90.]
        with pytest.raises(ValueError, match=r"All box lengths and angles"):
            box_matrix = convert_to_boxmatrix(box)

        box = [5., 5., 5., 182., 90., 90.]
        with pytest.raises(ValueError, match=r"All box angles"):
            box_matrix = convert_to_boxmatrix(box)

    def test_cubic(self):
        box = [5., 5., 5., 90., 90., 90.]
        box_matrix = convert_to_boxmatrix(box)
        assert box_matrix[0][0] == 5.0
        assert box_matrix[1][1] == 5.0
        assert box_matrix[2][2] == 5.0
        assert box_matrix[0][1] == 0.0
        assert box_matrix[0][2] == 0.0
        assert box_matrix[1][0] == 0.0
        assert box_matrix[1][2] == 0.0
        assert box_matrix[2][0] == 0.0
        assert box_matrix[2][1] == 0.0

    def test_rectangular(self):
        box = [5., 10., 15., 90., 90., 90.]
        box_matrix = convert_to_boxmatrix(box)
        assert box_matrix[0][0] == 5.0
        assert box_matrix[1][1] == 10.0
        assert box_matrix[2][2] == 15.0
        assert box_matrix[0][1] == 0.0
        assert box_matrix[0][2] == 0.0
        assert box_matrix[1][0] == 0.0
        assert box_matrix[1][2] == 0.0
        assert box_matrix[2][0] == 0.0
        assert box_matrix[2][1] == 0.0

    def test_triclinic(self):
        box = [10., 10., 10., 30., 75., 80.]
        box_matrix = convert_to_boxmatrix(box)
        assert box_matrix[0][0] == 10.0
        assert box_matrix[0][1] == 0.0
        assert box_matrix[0][2] == 0.0
        assert box_matrix[1][2] == 0.0
        assert np.isclose(box_matrix[1][1], 9.848078)
        assert np.isclose(box_matrix[2][2], 4.877255)
        assert np.isclose(box_matrix[1][0], 1.736482)
        assert np.isclose(box_matrix[2][0], 2.588191)
        assert np.isclose(box_matrix[2][1], 8.337484)

    def test_invalid_triclinic(self):
        box = [10., 10., 10., 60., 5., 90.]
        with pytest.raises(ValueError, match=r"Illegal box"):
            box_matrix = convert_to_boxmatrix(box)

