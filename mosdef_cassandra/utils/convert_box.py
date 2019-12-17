
import numpy as np

# NOTE: CODE COPIED FROM MDANALYSIS --> Need to fix

def triclinic_vectors(dimensions):
    """Convert ``[lx, ly, lz, alpha, beta, gamma]`` to a
    box matrix representation.

    Original `code by Tsjerk Wassenaar`_ posted on the Gromacs mailinglist.

    .. _code by Tsjerk Wassenaar:
       http://www.mail-archive.com/gmx-users@gromacs.org/msg28032.html

    Parameters
    ----------
    dimensions : list
        box dimensions
        ``[lx, ly, lz, alpha, beta, gamma]``.

    Returns
    -------
    box_matrix : numpy.ndarray
        A numpy array of shape ``(3, 3)`` with ``box_matrix[0]``
        containing the first, ``box_matrix[1]`` the second, and
        ``box_matrix[2]`` the third box vector.

    Notes
    -----
    * The first vector is guaranteed to point along the x-axis, i.e., it has the
      form ``(lx, 0, 0)``.
    * The second vector is guaranteed to lie in the x/y-plane, i.e., its
      z-component is guaranteed to be zero.
    * If any box length is negative or zero, or if any box angle is zero, the
      box is treated as invalid and an all-zero-matrix is returned.
    """

    dim = np.asarray(dimensions, dtype=np.float64)
    lx, ly, lz, alpha, beta, gamma = dim

    # Only positive edge lengths and angles in (0, 180) are allowed:
    if not np.all(dim > 0.0):
        raise ValueError('Invalid box dimensions. All box lengths'
                'must be > 0')
    if alpha > 180.0 or beta > 180.0 or gamma > 180.0:
        raise ValueError('Invalid box dimensions. All box angles'
                'must be > 0')

    # detect orthogonal boxes:
    if alpha == beta == gamma == 90.0:
        # box is orthogonal, return a diagonal matrix:
        box_matrix = np.diag(dim[:3])
    # we have a triclinic box:
    else:
        box_matrix = np.zeros((3, 3), dtype=np.float64)
        box_matrix[0, 0] = lx
        # Use exact trigonometric values for right angles:
        if alpha == 90.0:
            cos_alpha = 0.0
        else:
            cos_alpha = np.cos(np.deg2rad(alpha))
        if beta == 90.0:
            cos_beta = 0.0
        else:
            cos_beta = np.cos(np.deg2rad(beta))
        if gamma == 90.0:
            cos_gamma = 0.0
            sin_gamma = 1.0
        else:
            gamma = np.deg2rad(gamma)
            cos_gamma = np.cos(gamma)
            sin_gamma = np.sin(gamma)
        box_matrix[1, 0] = ly * cos_gamma
        box_matrix[1, 1] = ly * sin_gamma
        box_matrix[2, 0] = lz * cos_beta
        box_matrix[2, 1] = lz * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
        box_matrix[2, 2] = np.sqrt(lz * lz - box_matrix[2, 0] ** 2 - \
                                   box_matrix[2, 1] ** 2)
        # The discriminant of the above square root is only negative or zero for
        # triplets of box angles that lead to an invalid box (i.e., the sum of
        # any two angles is less than or equal to the third).
        # We don't need to explicitly test for np.nan here since checking for a
        # positive value already covers that.
        if box_matrix[2, 2] > 0.0:
            # all good, convert to correct dtype:
            box_matrix = box_matrix.astype(dtype=np.float64, copy=False)
        else:
            raise ValueError('Illegal box')
    return box_matrix
