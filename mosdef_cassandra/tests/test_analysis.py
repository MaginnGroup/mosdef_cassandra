import sys
import pytest
import numpy as np

from mosdef_cassandra.analysis import ThermoProps
from mosdef_cassandra.tests.base_test import BaseTest

try:
    import pandas as pd

    has_pandas = True
except ImportError:
    has_pandas = False


class TestAnalysis(BaseTest):
    def test_read_prp(self):
        thermo = ThermoProps("files/equil.out.box1.prp")
        assert thermo._data.shape == (201, 7)
        assert np.isclose(thermo._data[-1, 3], 42.242944)
        assert np.isclose(thermo._data[0, 2], 22378.33)

    def test_extract_prop(self):
        thermo = ThermoProps("files/equil.out.box1.prp")
        assert thermo.prop("Pressure").shape == (201,)
        assert np.isclose(thermo.prop("Pressure")[-1].value, 42.242944)
        assert np.isclose(thermo.prop("Energy_InterVDW")[0].value, 22378.33)

    def test_invalid_prop(self):
        thermo = ThermoProps("files/equil.out.box1.prp")
        with pytest.raises(ValueError, match=r"not an available"):
            chem_pot = thermo.prop("Chemical_Potential")

    def test_extract_range(self):
        thermo = ThermoProps("files/equil.out.box1.prp")
        pressure = thermo.prop("Pressure", start=15, end=30)
        sweep = thermo.prop("MC_SWEEP", start=15, end=30)

        assert pressure.shape == sweep.shape == (4,)
        assert np.isclose(sweep[0].value, 15)
        assert np.isclose(sweep[-1].value, 30)
        assert np.isclose(pressure[0].value, 152.82276)
        assert np.isclose(pressure[-1].value, -2.9842435)

    @pytest.mark.skipif(not has_pandas, reason="pandas not installed")
    def test_to_df(self):
        thermo = ThermoProps("files/equil.out.box1.prp")
        df = thermo.to_df()
        arrays = [
            (
                "MC_SWEEP",
                "Energy_Total",
                "Energy_InterVDW",
                "Pressure",
                "Volume",
                "Nmols",
                "Mass_Density",
            ),
            (
                "",
                "(kJ/mol)-Ext",
                "(kJ/mol)-Ext",
                "(bar)",
                "(A^3)",
                "",
                "(kg/m^3)",
            ),
        ]

        multi_index = pd.MultiIndex.from_arrays(
            arrays, names=("property", "units")
        )
        assert (df.columns == multi_index).all()
        assert df.shape == (201, 7)
