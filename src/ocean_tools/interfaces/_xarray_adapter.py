import xarray as xr
import dask.array as da
import numpy as np
from typing import Dict, Optional, Callable, List, Union, Any
from docstring_parser import parse, compose, DocstringParam, DocstringReturns
from datetime import datetime
import inspect

from ..version import __version__

def xarray_adapter(
    *input_fields_names,
    xarray=False,
    pass_self=False,
    expansion_dims=None,
    **output_fields_dims):
    """
    Decorate a function to work in a xarray Dataset.

    Given a function that works on dask or numpy arrays and generates one
    or multiple outputs, it is possible to specify that the inputs and outputs
    are respectively taken from and registered in the Dataset. This adapter
    takes the input fields names as *args and the outputs fields names and dimensions
    as **kwargs. This allows to generate a decorator that will work with an xarray Dataset.

    It is important that the input fields of the function can be passed either as a
    positional argument or a keyword argument (/ notation in function), because the decorator
    will systematically call the wrapped function using keyword arguments for the input fields.
    Moreover, the output fields of the wrapped function shall be returned as a dictionary to 
    allow identification of the output fields.


    Parameters
    ----------
    input_fields_names
        Input fields to extract from the dataset
    xarray
        Whether the wrapped function works with xarray or numpy arrays
    pass_self
        Whether the wrapped function is a class method and need self as a positional argument
    output_fields_dims
        Tuple with key as the output field name to register in the dataset and value as the
        dimensions in the dataset
    expansion_dims    
        Optional list of dimension over which to expand all variables. This is
        useful if the wrapped numpy function will use broadcasting
        functionalities between arrays of different shapes

    Returns
    -------
    :
        Decorator
    """

    def decorator(func):

        def wrapper(self, *args, **kwargs):

            # The user can override the input and output field named
            input_fields, output_fields = _override_fields_names(
                input_fields_names, output_fields_dims, kwargs
            )

            # Use the first input field to guess the function axes
            first_field_name = input_fields[input_fields_names[0]]
            kwargs.update(
                _populate_axis_identifiers(func, self._obj, first_field_name, self._axis_mapping, xarray=xarray)
            )

            # Initialize attributes.
            # This must be done before putting the dataset input arrays in the kwargs
            attributes = _initialize_attributes(
                input_fields, output_fields, func, args, kwargs
            )

            # Add fields in kwargs for function call
            kwargs.update(
                _extract_fields_from_dataset(
                    self._obj,
                    input_fields,
                    func,
                    xarray=xarray,
                    expansion_dims=expansion_dims))

            # Function call
            if pass_self:
                results = func(self, *args, **kwargs)
            else:
                results = func(*args, **kwargs)

            # Update underlying store
            self._obj = _update_dataset(self._obj, attributes, output_fields, results)
            return self


        wrapper.__doc__ = _edit_docstring(func, input_fields_names, output_fields_dims.keys())
        return wrapper

    return decorator

def _update_dataset(
    dataset: xr.Dataset,
    attributes: Dict[str, str],
    output_fields: Dict[str, Any],
    results: Dict[str, np.ndarray],
) -> xr.Dataset:

    # Convert to dict if there is only one output
    if not isinstance(results, dict):
        result = results
        results = dict()
        results[list(output_fields.keys())[0]]=result

    # Create output xarray.DataArrays, with attributes to track the processing
    arrays=dict()
    for k, v in results.items():
        arrays[output_fields[k]["new_name"]] = xr.DataArray(
            data=v,
            attrs=attributes[k],
            dims=output_fields[k]["dims"]
        )

    return dataset.update(arrays)

def _extract_fields_from_dataset(
    dataset: xr.Dataset,
    input_fields: Dict[str, str],
    func: Callable,
    xarray: bool = False,
    expansion_dims: Optional[List[str]] = None,
) -> Union[Dict[str, np.ndarray], Dict[str, xr.DataArray]]:
    extracted_fields = dict()
    for k, v in input_fields.items():
        # Retrieve field in dataset and rewrite the error message for clarity
        try:
            extracted_field = dataset[v]
        except:
            raise Exception(
                f"Function {func.__name__} needs input_field '{v}'"
                f" but it does not exist in dataset").with_traceback(None) from None

        # Expand field if necessary
        if expansion_dims is not None:
            missing_dims = { d: ii for ii, d in enumerate(expansion_dims) if d not in extracted_field.dims}
            
            extracted_field = extracted_field.expand_dims(
                list(missing_dims.keys()), list(missing_dims.values()))

        # Get either numpy or xarray array
        if not xarray:
            extracted_field = extracted_field.data

        # Put in dictionary        
        extracted_fields[k] = extracted_field
    
    return extracted_fields

def _populate_axis_identifiers(
    func: Callable,
    dataset: xr.Dataset,
    first_field_name: str,
    axis_mapping: Dict[str, int],
    xarray: bool = False,
) -> Union[Dict[str, str], Dict[str, int]]:
    parameters = inspect.signature(func).parameters.keys()
    try:
        first_field = dataset[first_field_name]
    except KeyError:
        raise Exception(
            f"Function {func.__name__} needs input_field '{first_field_name}'"
            " but it does not exist in dataset").with_traceback(None) from None

    # An axis can have multiple dimensions associated. Exemple MITGCM has an X axis with the
    # centered ticks (i) and the left ticks (i_g)
    unique_axis_mapping = dict()
    for k, v in axis_mapping.items():
        if isinstance(v, list):
            unique_axis_mapping[k] = [vv for vv in v if vv in first_field.dims][0]
        else:
            unique_axis_mapping[k] = v

    axes = dict()
    for k, v in unique_axis_mapping.items():
        if k in parameters:
            try:
                axes[k] = v if xarray else first_field.dims.index(v)
            except ValueError:
                raise Exception(
                    f"Function {func.__name__} needs to index axis '{v}' but it does not exist"
                    f" in the first variable '{first_field.name}' used to extract the index (dims={first_field.dims})")

    return axes


def _override_fields_names(
    input_fields_names,
    output_fields_dims,
    kwargs
):
    input_fields = {name: kwargs.pop(name) if name in kwargs else name for name in input_fields_names}
    output_fields = dict()
    for k, v in output_fields_dims.items():
        output_fields[k] = dict(
            new_name=kwargs.pop(k) if k in kwargs else k,
            dims=v,
        )

    return input_fields, output_fields


def _initialize_attributes(
    input_fields,
    output_fields,
    process_func,
    process_args,
    process_kwargs,
):
    # Processing information
    attributes_template = dict()
    attributes_template["function"] = process_func.__qualname__
    attributes_template["inputs_fields"] = str(list(input_fields.values()))
    if len(process_args) > 0:
        attributes_template["args"] = str(process_args)
    if len(process_kwargs) > 0:
        attributes_template["kwargs"] = str(process_kwargs)
    attributes_template["date"] = datetime.now().strftime("%Y-%m-%dT%H:%M:%S.%f")
    attributes_template["sofware_version"] = __version__

    # Parse current docstring. This gives us access to the fields descriptions
    doc_string = parse(process_func.__doc__)
    doc_string_returns = {r.return_name: r.description for r in doc_string.many_returns}

    attributes = dict()
    for k in output_fields.keys():
        try:
            description = doc_string_returns[k]
        except KeyError:
            print(f"Warning: could not find description of return fields {k} in function docstring")
            description = None

        attributes[k] = attributes_template.copy()
        attributes[k]["description"] = description

    return attributes


def _edit_docstring(func: Callable, input_fields: List[str], output_fields: List[str]):
    """
    Change the docstring of a function wrapped with xarray_adapter decorator.
    
    A function wrapped to work with an xarray dataset will extract some of its inputs
    from a dataset, and will write some of its outputs in the same dataset. As such,
    inputs that were arrays are now names in the dataset. The same goes for the
    outputs, but in addition their name are now inputs. All of this must be reflected
    in the docstring.
    
    Parameters
    ----------
    input_fields
        List of input fields extracted from the dataset
    output_fields
        List of output fields registered in the dataset
    
    Returns
    -------
    :
        Edited docstring
    """
    doc_string = parse(func.__doc__)
    doc_string_returns = {r.return_name: r.description for r in doc_string.many_returns}

    # Update function docstring
    for param in doc_string.params:
        if param.arg_name in input_fields:
            param.description += " (Input name in dataset)"
            param.is_optional = True
            param.type_name = "str"

    for k in output_fields:
        try:
            doc_string.meta.append(
                DocstringParam(["param", k], doc_string_returns[k] + " (Output name in dataset)", k, "str", True, k)
            )
        except:
            print(f"Warning: field {k} does not exist in the original docstring of {func.__qualname__}. Docstring might be malformed")

    # Remove returned values (they are set in the dataset)
    doc_string.meta = [meta for meta in doc_string.meta if not isinstance(meta, DocstringReturns)]
    doc_string.meta.append(
        DocstringReturns(
            ["returns", "xarray_accessor"],
            "Original instance of the xarray accessor, which contains the xarray dataset populated with the results",
            "Accessor",
            False
        )
    )
        
    return compose(doc_string)