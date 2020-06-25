import pytest
import unyt as u
import mosdef_cassandra as mc

from unyt import dimensions
from unyt.array import allclose_units
from mosdef_cassandra.tests.base_test import BaseTest
from mosdef_cassandra.writers.inp_functions import (
    _check_kwarg_units_helper,
    _convert_kwarg_units_helper,
)
from mosdef_cassandra.utils.units import validate_unit_list


class TestKwargUnits(BaseTest):
    def test_check_kwarg_units_helper(self):
        kwargs = {
            "example_cutoff": 5.5 * u.angstrom,
            "example_cutoff_list": [5.5 * u.angstrom, 6.0 * u.angstrom],
        }
        _check_kwarg_units_helper(kwargs, "example_cutoff", dimensions.length)
        _check_kwarg_units_helper(
            kwargs, "example_cutoff", dimensions.length, list_length=2
        )
        with pytest.raises(TypeError, match=""):
            _check_kwarg_units_helper(
                kwargs, "example_cutoff", dimensions.pressure
            )
        _check_kwarg_units_helper(
            kwargs, "example_cutoff_list", dimensions.length, list_length=2
        )
        with pytest.raises(TypeError, match=""):
            _check_kwarg_units_helper(
                kwargs,
                "example_cutoff_list",
                dimensions.pressure,
                list_length=2,
            )
        with pytest.raises(TypeError, match=""):
            _check_kwarg_units_helper(
                kwargs, "example_cutoff_list", dimensions.length, list_length=3
            )
        with pytest.raises(TypeError, match=""):
            _check_kwarg_units_helper(
                kwargs, "example_cutoff_list", dimensions.length
            )

    def test_convert_kwarg_units_helper(self):
        kwargs = {
            "example_cutoff": 5.5 * u.angstrom,
            "example_cutoff_list": [5.5 * u.angstrom, 6.0 * u.angstrom],
            "example_cutoff_array": [5.5, 6.0] * u.angstrom,
        }
        _convert_kwarg_units_helper(kwargs, "example_cutoff", "nm")
        assert kwargs["example_cutoff"].units == u.nm
        assert allclose_units(kwargs["example_cutoff"], 0.55 * u.nm)
        with pytest.raises(TypeError, match="went wrong"):
            _convert_kwarg_units_helper(kwargs, "example_cutoff_list", "nm")
        _convert_kwarg_units_helper(kwargs, "example_cutoff_array", "nm")
        assert kwargs["example_cutoff_array"].units == u.nm
        assert allclose_units(
            kwargs["example_cutoff_array"], [0.55, 0.60] * u.nm
        )

    def test_validate_unit_list(self):
        list_ = [3.0 * u.angstrom, 1.0 * u.angstrom, 1.5 * u.angstrom]
        with pytest.raises(TypeError):
            validate_unit_list(list_, (2,), dimensions.length)
