from datetime import datetime
import functools
import inspect
#from ocean_tools.version import __version__

def stamp_processing_info(xarray_func):
    """Add processing date and software version used to edit some datasets fields.
    
    This wrapper only works on a function that returns a tuple of xarray
    DataArrays.
    """
    @functools.wraps(xarray_func)
    def wrapped(*args, **kwargs):

        arrays = xarray_func(*args, **kwargs)

        if not isinstance(arrays, tuple):
            single_output = True
            arrays = (arrays,)
        else:
            single_output=False
        
        new_arrays = []
        for array in arrays:
            array.attrs["processing"] = _processing_info(xarray_func, *args, **kwargs)
            new_arrays.append(array)
        
        if single_output:
            return new_arrays[0]
        else:
            return tuple(new_arrays)

    return wrapped


def _processing_info(func, *args, **kwargs):
    
    processing = dict()
    processing["last_modified_date"] = datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")
    #processing["sofware_version"] = __version__
    processing["function"] = func.__qualname__

    return processing


def _input_attributes(func, *args, **kwargs):
    # Get default arguments
    # TODO: this will help the processing traceability, but we need a more
    # complex system than just putting that in the attributes
    arguments = dict()
    for p in inspect.signature(func).parameters.values():
        if p.default is not p.empty:
            arguments[p.name] = p.default
    arguments.update(kwargs)

    # Extract positional arguments of a function
    inputs = {
        k.name if hasattr(k, "name") and k.name is not None else f"arg{ii}": k.attrs if hasattr(k, "attrs") else "<Positional argument without attributes>"
        for ii, k in enumerate(args)}

    # Extract keyword arguments
    inputs.update(
        {k: v.attrs if hasattr(v, "attrs") else str(v)
        for k, v in arguments.items()})

    return inputs
