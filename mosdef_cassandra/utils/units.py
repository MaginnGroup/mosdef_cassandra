import unyt as u
from unyt import accepts


def validate_unit(unyt_array, valid_dimension, argument_name):
    if isinstance(unyt_array, u.unyt_array):

        @accepts(unyt_array=valid_dimension)
        def _validate(unyt_array):
            return unyt_array

    else:
        raise TypeError(f"{argument_name} must be a unyt_array")

    try:
        return _validate(unyt_array)
    except TypeError:
        raise TypeError(
            f"{argument_name} must be a list of unyt_quantities or "
            f"unyt_array with dimensions of {valid_dimension}"
        )


def validate_unit_list(list_, valid_shape, valid_dimension, argument_name):
    try:
        unyt_array = u.unyt_array(list_)
    except AttributeError:
        raise TypeError(
            f"{argument_name} must be a list of unyt_quantities or "
            f"unyt_array with shape = {valid_shape} and dimensions "
            f"of {valid_dimension}"
        )
    if unyt_array.shape != valid_shape:
        raise TypeError(
            f"{argument_name} must be a list of unyt_quantities or "
            f"unyt_array with shape = {valid_shape} and dimensions "
            f"of {valid_dimension}"
        )

    return validate_unit(unyt_array, valid_dimension, argument_name)
