import pytest
import numpy as np
import unyt as u
import mbuild

from mosdef_cassandra.tests.base_test import BaseTest
from mosdef_cassandra.utils.units import validate_unit, validate_unit_list
from unyt import dimensions
from unyt.exceptions import IterableUnitCoercionError


class TestConvertBox(BaseTest):
    def test_cubic(self):
        box = mbuild.Box([5.0, 5.0, 5.0], [90.0, 90.0, 90.0])
        box_matrix = box.vectors
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
        box = mbuild.Box([5.0, 10.0, 15.0], [90.0, 90.0, 90.0])
        box_matrix = box.vectors
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
        box = mbuild.Box([10.0, 10.0, 10.0], [30.0, 75.0, 80.0])
        box_matrix = box.vectors
        assert box_matrix[0][0] == 10.0
        assert box_matrix[0][1] == 0.0
        assert box_matrix[0][2] == 0.0
        assert box_matrix[1][2] == 0.0
        assert np.isclose(box_matrix[1][1], 9.848078)
        assert np.isclose(box_matrix[2][2], 4.877255)
        assert np.isclose(box_matrix[1][0], 1.736482)
        assert np.isclose(box_matrix[2][0], 2.588191)
        assert np.isclose(box_matrix[2][1], 8.337484)

    @pytest.mark.parametrize(
        "unit,dimension,name",
        [
            (u.nm, dimensions.length, "length"),
            (u.bar, dimensions.pressure, "pressure"),
            (u.K, dimensions.temperature, "temperature"),
            ((u.kJ / u.mol), dimensions.energy, "energy"),
        ],
    )
    def test_validate_unit(self, unit, dimension, name):
        validate_unit(1 * unit, dimension, argument_name=name)

    def test_validate_unit_int_error(self):
        with pytest.raises(TypeError):
            validate_unit(1, dimensions.length)

    def test_invalid_dimension(self):
        with pytest.raises(TypeError, match="with dimensions of"):
            validate_unit(1 * u.nm, dimensions.temperature)

    def test_unit_err_msg(self):
        with pytest.raises(TypeError, match="test must be a"):
            validate_unit(
                1 * u.nm, dimensions.temperature, argument_name="test"
            )

    @pytest.mark.parametrize(
        "unit_list, shape, dimension",
        [
            ([1.0 * u.nm], (1,), dimensions.length),
            ([[1.0 * u.nm]], (1, 1), dimensions.length),
            ([[1.0 * u.nm, 1.0 * u.nm]], (1, 2), dimensions.length),
            ([[1.0, 1.0] * u.nm], (1, 2), dimensions.length),
            ([[1.0, 1.0]] * u.nm, (1, 2), dimensions.length),
            ([[1.0 * u.nm], [1.0 * u.nm]], (2, 1), dimensions.length),
            ([[1.0] * u.nm, [1.0] * u.nm], (2, 1), dimensions.length),
            ([[1.0], [1.0]] * u.nm, (2, 1), dimensions.length),
            (
                [[1.0, 1.0] * u.nm, [1.0 * u.nm, 1.0 * u.nm]],
                (2, 2),
                dimensions.length,
            ),
        ],
    )
    def test_validate_unit_list(self, unit_list, shape, dimension):
        unit_list = validate_unit_list(unit_list, shape, dimension)
        assert type(unit_list) == u.unyt_array
        assert unit_list.units.dimensions == dimension
        assert unit_list.shape == shape

    @pytest.mark.parametrize(
        "unit_list, shape, dimension",
        [
            ([1.0 * u.nm], (), dimensions.length),
            ([[1.0 * u.nm]], (2, 1), dimensions.length),
            ([[1.0 * u.nm, 1.0 * u.nm]], (2, 1), dimensions.length),
            ([[1.0, 1.0] * u.nm], (1, 2), dimensions.pressure),
            ([[1.0, 1.0 * u.nm]], (2, 1), dimensions.length),
        ],
    )
    def test_invalid_unit_list(self, unit_list, shape, dimension):
        with pytest.raises(TypeError, match="argument must be a list"):
            validate_unit_list(unit_list, shape, dimension)
