import numpy as np


def convert_to_boxmatrix(dimensions):
    """Convert ``[x, y, z, alpha, beta, gamma]`` to a
    matrix representation.

    Original code by Tsjerk Wassenaar posted on Gromacs mailinglist.

       http://www.mail-archive.com/gmx-users@gromacs.org/msg28032.html

    Parameters
    ----------
    dimensions : list
        box dimensions
        ``[x, y, z, alpha, beta, gamma]``.

    Returns
    -------
    box_matrix : numpy.ndarray
        A numpy array of shape ``(3, 3)`` with ``box_matrix[0]``
        containing the first, ``box_matrix[1]`` the second, and
        ``box_matrix[2]`` the third box vector.

    Notes
    -----
    * The first box vector points along the x-axis, i.e., it has the
      form ``(x, 0, 0)``.
    * The second box vector lies in the x/y-plane, i.e., its has the
      form ``(x', y', 0)``.
    * All box lengths must be positive.
    """

    dim = np.asarray(dimensions)
    # Sanity checks
    if dim.shape != (6,):
        raise ValueError(
            "Invalid box dimensions. Input must be provided as "
            "[x, y, z, alpha, beta, gamma]"
        )
    if not np.all(dim > 0.0):
        raise ValueError(
            "Invalid box dimensions. All box lengths and angles must be > 0"
        )

    x, y, z, alpha, beta, gamma = dim

    # Only angles in (0, 180) are allowed:
    if alpha > 180.0 or beta > 180.0 or gamma > 180.0:
        raise ValueError(
            "Invalid box dimensions. All box angles must be < 180 degrees"
        )

    box_matrix = np.zeros((3, 3))
    # First check for orthogonal boxes
    if alpha == beta == gamma == 90.0:
        box_matrix[0][0] = x
        box_matrix[1][1] = y
        box_matrix[2][2] = z
    # Else triclinic
    else:
        # Calc trig values
        cos_alpha = np.cos(np.radians(alpha))
        cos_beta = np.cos(np.radians(beta))
        cos_gamma = np.cos(np.radians(gamma))
        sin_gamma = np.sin(np.radians(gamma))
        # Protect against numerical errors for 90 deg
        cos_alpha = np.round(cos_alpha, 15)
        cos_beta = np.round(cos_beta, 15)
        cos_gamma = np.round(cos_gamma, 15)
        sin_gamma = np.round(sin_gamma, 15)
        # Calc box dims
        box_matrix[0][0] = x
        box_matrix[1][0] = y * cos_gamma
        box_matrix[1][1] = y * sin_gamma
        box_matrix[2][0] = z * cos_beta
        box_matrix[2][1] = z * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
        box_matrix[2][2] = np.sqrt(
            z * z - box_matrix[2, 0] ** 2 - box_matrix[2, 1] ** 2
        )
        if np.isnan(box_matrix[2][2]):
            raise ValueError("Illegal box")
    return box_matrix
