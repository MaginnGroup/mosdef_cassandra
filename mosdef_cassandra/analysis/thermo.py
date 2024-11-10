import os
import numpy as np
import unyt as u


class ThermoProps:
    """Store thermodynamic properties from a Cassandra .prp file"""

    def __init__(self, filename):
        """Create ThermoProps from a .prp file

        Parameters
        ----------
        filename : string
            path to the .prp file

        Returns
        -------
        ThermoProps
            object containing the contents of the .prp file
        """

        self.filename = filename

        # Read headers
        prp_headers = []
        with open(self.filename) as f:
            for idx, line in enumerate(f):
                if idx == 1 or idx == 2:
                    prp_headers.append(line)
                if idx > 2:
                    break

        # Extract column names and units
        column_names = prp_headers[0][1:].split()
        column_units = [""]
        n_columns = len(column_names)
        for icol in range(1, n_columns):
            col_start = 12 + 18 * (icol - 1)
            col_end = 12 + 18 * icol
            column_units.append(prp_headers[1][col_start:col_end].strip())

        prp_data = np.genfromtxt(self.filename, skip_header=3)

        assert prp_data.shape[1] == len(column_names)
        assert prp_data.shape[1] == len(column_units)

        # Unyt mapping

        units_to_unyts = {
            "(kJ/mol)-Ext": u.Unit("kJ/mol"),
            "(bar)": u.bar,
            "": u.dimensionless,
            "(A^3)": u.angstrom**3,
            "(kg/m^3)": u.Unit("kg/m**3"),
            "(molec/A^3)": u.Unit("count/angstrom**3"),
        }

        self._properties = column_names
        unyts = []
        for unit in column_units:
            try:
                unyts.append(units_to_unyts[unit])
            except KeyError:
                unyts.append(u.dimensionless)

        self._units = column_units
        self._unyts = unyts
        self._data = prp_data

    def print_props(self):
        """Print the available properties"""
        for prop in self._properties:
            print(prop)

    def prop(self, prp_name, start=None, end=None, units=True):
        """Extract the specified property

        Parameters
        ----------
        prp_name : string
            the property to extract
        start : int
            the starting step/sweep/etc.
        end : int
            the ending step/sweep/etc.

        Returns
        -------
        unyt_array
            the property with units
        """
        if prp_name not in self._properties:
            raise ValueError(
                "{} is not an available property. Please select from: {}".format(
                    prp_name, self._properties
                )
            )
        col_idx = self._properties.index(prp_name)

        if start is not None:
            start_idx = np.where(self._data[:, 0] >= start)[0][0]
        else:
            start_idx = None
        if end is not None:
            end_idx = np.where(self._data[:, 0] <= end)[0][-1] + 1
        else:
            end_idx = None

        return self._data[start_idx:end_idx, col_idx] * self._unyts[col_idx]

    def to_df(self):
        """Convert ThermoProps to a pandas.DataFrame"""
        try:
            import pandas as pd
        except ModuleNotFoundError:
            raise ModuleNotFoundError(
                "The pandas package is required to convert to a pandas.DataFrame. "
                "pandas can be installed with 'conda install -c conda-forge pandas'"
            )

        arrays = [(self._properties), (self._units)]
        multi_index = pd.MultiIndex.from_arrays(
            arrays, names=("property", "units")
        )

        return pd.DataFrame(self._data, columns=multi_index)

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, filename):
        if os.path.exists(filename):
            self._filename = filename
        else:
            raise FileNotFoundError(
                "File {} could not be found".format(filename)
            )
