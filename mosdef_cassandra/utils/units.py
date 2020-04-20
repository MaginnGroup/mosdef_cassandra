import unyt as u
from unyt import accepts

def validate_unit(unyt_array, dimension):
    @accepts(unyt_array=dimension)
    def validate_unit(unyt_array):
        if isinstance(unyt_array, u.unyt_array):
            return unyt_array
        else:
            raise TypeError('Argument must be a unyt_array')

    return validate_unit(unyt_array)
