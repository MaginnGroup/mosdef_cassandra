import unyt as u
from unyt import accepts


def validate_unit(
    unyt_array, valid_dimension, argument_name=None, err_msg=None
):

    if argument_name is None:
        argument_name = "argument"
    if err_msg is None:
        err_msg = (
            f"{argument_name} must be a unyt_quantity or unyt_array "
            f"with dimensions of {valid_dimension}"
        )

    if isinstance(unyt_array, u.unyt_array):

        @accepts(unyt_array=valid_dimension)
        def _validate(unyt_array):
            return unyt_array

    else:
        raise TypeError(err_msg)

    try:
        return _validate(unyt_array)
    except TypeError:
        raise TypeError(err_msg)


def validate_unit_list(
    list_, valid_shape, valid_dimension, argument_name=None, err_msg=None
):

    if argument_name is None:
        argument_name = "argument"
    if err_msg is None:
        err_msg = (
            f"{argument_name} must be a list of unyt_quantities or "
            f"unyt_array with shape = {valid_shape} and dimensions "
            f"of {valid_dimension}"
        )

    if type(list_) == list:
        new_list = []
        for item in list_:
            try:
                shape = (len(item),)
            except TypeError:
                shape = ()
            new_list.append(
                validate_unit_list(
                    item, shape, valid_dimension, argument_name, err_msg
                )
            )
        list_ = new_list
    try:
        unyt_array = u.unyt_array(list_)
    except AttributeError:
        raise TypeError(err_msg)
    if unyt_array.shape != valid_shape:
        raise TypeError(err_msg)

    return validate_unit(unyt_array, valid_dimension, argument_name, err_msg)
