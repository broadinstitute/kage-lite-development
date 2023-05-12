import logging
import types

import numpy as np


def replace_object_attributes_recursively(object, func, ignore_types=None, key=None):
    """
    Runs replace_func on an object and its "children" recursively.
    Can be used to e.g. search for specific data and replace it.
    """

    if ignore_types is not None and isinstance(object, ignore_types):
        # dont unpack types in ignore_types
        return func(object, key=key)
    elif object is None:
        return object
    elif isinstance(object, (str, int, float, np.integer)):
        # don't do anything with primitive types
        return func(object, key=key)
    elif isinstance(object, np.ndarray):
        # don't do anything with numpy arrays
        res = func(object, key=key)
        return res
    elif isinstance(object, tuple):
        # convert to list, replace, and convert back
        return tuple(replace_object_attributes_recursively(list(object), func, ignore_types, key=key))
    elif isinstance(object, set):
        # just pickle set
        return object
    elif isinstance(object, dict):
        # run on each value
        return {k: replace_object_attributes_recursively(value, func, ignore_types, key=key + '_' + k if key is not None else k) for k, value in object.items()}
    elif isinstance(object, list):
        # run on each element
        return [replace_object_attributes_recursively(e, func, ignore_types, key=key + '_' + str(i) if key is not None else str(i)) for i, e in enumerate(object)]
    elif isinstance(object, np.dtype):
        return object
    elif callable(object):   #isinstance(object, types.FunctionType):
        #print("Object is callable: ", object)
        return object
    else:
        # assume regular object, iterate all attributes
        try:
            items = object.__dict__.items()
        except AttributeError:
            print(object)
            raise

        for attr, value in items:
            setattr(object, attr, replace_object_attributes_recursively(value, func, ignore_types, key=key + '_' + attr if key is not None else attr))

        return func(object, key=key)

