import pytest
import numpy as np
import unyt as u

from mosdef_cassandra.tests.base_test import BaseTest
from mosdef_cassandra.utils.convert_box import convert_to_boxmatrix
from mosdef_cassandra.utils.units import validate_unit
from unyt import dimensions


class TestConvertBox(BaseTest):
    def test_invalid_shape(self):
        box = [10.0, 5.0, 5.0]
        with pytest.raises(ValueError, match=r"Input must be provided as"):
            box_matrix = convert_to_boxmatrix(box)

    def test_invalid_lengths(self):
        box = [10.0, 5.0, -5.0, 90.0, 90.0, 90.0]
        with pytest.raises(ValueError, match=r"All box lengths and angles"):
            box_matrix = convert_to_boxmatrix(box)

    def test_invalid_angles(self):
        box = [5.0, 5.0, 5.0, -20.0, 90.0, 90.0]
        with pytest.raises(ValueError, match=r"All box lengths and angles"):
            box_matrix = convert_to_boxmatrix(box)

        box = [5.0, 5.0, 5.0, 182.0, 90.0, 90.0]
        with pytest.raises(ValueError, match=r"All box angles"):
            box_matrix = convert_to_boxmatrix(box)

    def test_cubic(self):
        box = [5.0, 5.0, 5.0, 90.0, 90.0, 90.0]
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
        box = [5.0, 10.0, 15.0, 90.0, 90.0, 90.0]
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
        box = [10.0, 10.0, 10.0, 30.0, 75.0, 80.0]
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
        box = [10.0, 10.0, 10.0, 60.0, 5.0, 90.0]
        with pytest.raises(ValueError, match=r"Illegal box"):
            box_matrix = convert_to_boxmatrix(box)

    @pytest.mark.parametrize(
        "unit,dimension",
        [
            (u.nm, dimensions.length),
            (u.bar, dimensions.pressure),
            (u.K, dimensions.temperature),
            ((u.kJ / u.mol), dimensions.energy),
        ],
    )
    def test_validate_unit(self, unit, dimension):
        validate_unit(1 * unit, dimension)

    def test_validate_unit_int_error(self):
        with pytest.raises(TypeError):
            validate_unit(1, dimensions.length)

    def test_invalid_dimension(self):
        with pytest.raises(TypeError, match="does not match"):
            validate_unit(1 * u.nm, dimensions.temperature)
