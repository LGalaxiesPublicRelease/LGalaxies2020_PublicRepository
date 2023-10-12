"""
Script containing functions for ensembles
"""

import operator
from functools import reduce  # forward compatibility for Python 3

import collections
from collections import OrderedDict

import numpy as np
import astropy.units as u

# from binarycpython.utils.dicts import AutoVivificationDict

# function to provide a list that can act as a sequence of key selections
def getFromDict(dataDict, mapList):
    return reduce(operator.getitem, mapList, dataDict)

def more_or_less(value, epsilon=1e-5):
    """
    Function to create a lower and upper bound around a given value
    """

    if isinstance(value, (int, float, np.float64, np.int64)):
        return [value-epsilon, value+epsilon]
    elif isinstance(value, (list, tuple, np.array)):
        return [value[0]-epsilon, value[-1]+epsilon]
    else:
        raise ValueError("Object type {} not supported.".format(type(value)))

def query_dict(ensemble_dict, keyname, valuerange, depth=0, max_depth=-1, use_max_depth=False, verbose=0):
    """
    Function to query a nested dictionary.

    This is a recursive function that will try to match the keyname.

    If it finds a match, it will try to match the keys whether they fall within the given valuerange.
        When that given key is found, it will add that subdict to the dict that will be returned, but it wont call itself anymore.

    If it does not match, it will recurse downward until the value is not a dict anymore.

    During the recursive walk we track the depth.
        If the function is configured to use_max_depth we will track which depth matches the key we want to match. TODO: consider if this is actually useful

    The best way to use this function is we know that the (depth)structure is the same. I.e. that on every depth the layer contains the same keys

    Input:
        ensemble_dict: the ensemble dict that we want to query
        keyname: key that we want to query.
        valuerange: range in which the value should fall. The value will be evaluated as a float(value) since they are generally stored as a string.
    Output:
        dict with matching subdicts, or an empty dict.
    """

    return_dict = {}

    for key in ensemble_dict:
        # if the key matches the desired keyname, then check its subdict values whether they match the values that we want
        if key == keyname:
            if verbose:
                print('matched key: {}'.format(key))
            return_subdict = {} # create empty subdict
            for subkey in ensemble_dict[key]:
                if valuerange[0] <= float(subkey) <= valuerange[1]:
                    # Fill subdict if its matched
                    if verbose:
                        print("subkey {} falls within range [{}, {}]".format(subkey, valuerange[0], valuerange[1]))
                    return_subdict[subkey] = ensemble_dict[key][subkey]

        # if not, check if the subdict is a dictionary type
        else:
            if not isinstance(ensemble_dict[key], (dict, OrderedDict)):
                return_subdict = {}
            # If it is, recurse
            else:
                return_subdict = query_dict(ensemble_dict[key], keyname, valuerange, depth=depth+1, max_depth=max_depth, use_max_depth=use_max_depth, verbose=verbose)

        # If the subdict that is found has some content, add it to the return_dict
        if return_subdict:
            return_dict[key] = return_subdict

    #
    return return_dict


def return_keys_bins_and_bincenters(input_ensemble, key):
    """
    Function to find all the keys belonging to the input_ensemble[key] subdict

    Then uses those to generate the bins from those keys (finding smallest difference and generating range of bins)

    Then also returns the centers of these bins
    """

    key_list = return_float_list_of_dictkeys(input_ensemble[key])
    bins = generate_bins_from_keys(key_list)
    bincenters = (bins[1:] + bins[:-1])/2

    return key_list, bins, bincenters

def return_float_list_of_dictkeys(input_dict):
    """
    function to return the sorted list where they keys are made floats
    """

    return np.array(sorted([float(el) for el in input_dict.keys()]))

def generate_bins_from_keys(input_array):
    """
    Function to generate the bins that the array of keys represents
    """

    stepsize = np.around(min(input_array[1:]-input_array[:-1]), decimals=10)
    bins = np.arange(
        input_array.min() - 0.5*stepsize,
        input_array.max() + 1.5*stepsize,
        stepsize
    )

    return bins


def flatten_data_ensemble2d(input_dict, named_subkey_1, named_subkey_2):
    """
    Function that loops over the nested dictionaries and creates a list where the first numerical subkey,
    the second numerical subkey and the actual value are combined in a list and then appended to a grand list.

    This requires the input dictionary to be structured as follows:

    {
        named_subkey_1:
        {
            numerical subkeys 1:
            {
                named_subkey_2:
                {
                    numerical_subkeys 2:
                        value,
                    ..
                }
            }
        }
    }

    TODO: turn this into a general function where the named subkeys are provided in a list that a recursive function can use
    """

    data = []
    for numerical_key_1 in sorted(input_dict[named_subkey_1]):
        for numerical_key_2 in input_dict[named_subkey_1][numerical_key_1][named_subkey_2]:
            value = input_dict[named_subkey_1][numerical_key_1][named_subkey_2][numerical_key_2]
            data.append([float(numerical_key_1), float(numerical_key_2), float(value)])
    data = np.array(data).T

    return data


def flatten_data_ensemble1d(input_dict, named_subkey_1):
    """
    Functon to get the subkey and its associated value from a dictionary.

    This dictionary can't be deeper than 2 levels (i.e. it has to be the last named key.

    This requires the input dictionary to be structured as follows:

    {
        named_subkey_1:
        {
            numerical subkeys 1:
                value,
            ..
        }
    }
    """

    data = []
    for numerical_key_1 in sorted(input_dict[named_subkey_1]):
        value = input_dict[named_subkey_1][numerical_key_1]

        try:
            converted_key = float(numerical_key_1)
        except:
            converted_key = str(numerical_key_1)
        data.append([converted_key, float(value)])
    data = np.array(data).T

    return data

#################
# Recursive counts and recursive sums: These recursive functions will go until the final depth and count or add the result
def get_recursive_count(input_dict):
    """
    Function to count the total amount of keys in a dictionary

    The ensemble dictionary can only contain dicts, lists or numbers.

    Numbers will get added, dicts will be called upon and lists will get skipped
    """

    local_count = 0
    if isinstance(input_dict, (dict, OrderedDict)):

        for key in input_dict.keys():
            local_count += 1
            if isinstance(input_dict[key], (dict, OrderedDict)):
                local_count += get_recursive_sum(input_dict[key])

            elif not isinstance(input_dict[key], list):
                local_count += input_dict[key]

    elif isinstance(input_dict, (int, float)):
        return input_dict

    else:
        raise ValueError("Unsupported type: {}".format(type(input_dict)))

    return local_count


def get_recursive_sum(input_dict):
    """
    Function to count the total amount of keys in a dictionary

    The ensemble dictionary can only contain dicts, lists or numbers.

    Numbers will get added, dicts will be called upon and lists will get skipped
    """

    local_count = 0
    if isinstance(input_dict, (dict, OrderedDict)):
        for key in input_dict.keys():
            if isinstance(input_dict[key], (dict, OrderedDict)):
                local_count += get_recursive_sum(input_dict[key])

            elif not isinstance(input_dict[key], list):
                local_count += input_dict[key]
    elif isinstance(input_dict, (int, float)):
        return input_dict

    else:
        raise ValueError("Unsupported type: {}".format(type(input_dict)))

    return local_count


def get_recursive_sum_for_subkeys_1d(input_dict, main_key):
    """
    Function that takes an input_dict[main_key] and loops over the keys.

    for the resulting input_dict[main_key][key] we call get_recursive_sum

    TODO: place example
    """

    new_dict = {}

    # Loop over keys:
    for key in input_dict[main_key]:
        new_dict[key] = get_recursive_sum(input_dict[main_key][key])

    return {main_key: new_dict}


def get_recursive_sum_for_subkeys_2d(input_dict, main_key, second_key):
    """
    Function that takes an input_dict[main_key] and loops over the keys.

    for the resulting input_dict[main_key][key][second_key] we call get_recursive_sum

    TODO: Turn this function into a more general function where the keys are a list that gets popped by a recursive function
    """

    new_dict = {}

    # Loop over keys:
    for key in input_dict[main_key]:

        new_dict[key] = get_recursive_sum_for_subkeys_1d(input_dict[main_key][key], second_key)        

    return {main_key: new_dict}

#################
# Merge all subdicts functions: Functions to (recursively) merge the subdicts

def merge_recursive_until_stop(ensemble_dict, stopping_key):
    """
    Function to merge the dictionary until we hit the stopping key
    """

    next_keys = list(ensemble_dict.keys())
    if not len(next_keys) == 1:
        raise ValueError("Things are going wrong. next_key: {}".format(next_keys))

    # Check if we arrived at the key
    next_key = next_keys[0]
    if next_key == stopping_key:
        return ensemble_dict

    # If we did not merge merge mergeee
    return merge_recursive_until_stop(merge_all_subdicts(ensemble_dict[next_key]), stopping_key)

def recursively_merge_all_subdicts(input_ensemble, key_list):
    """
    Function to recursively merge all the dicts:

    when calling this function the first time with:
        'recurse_call(input_ensemble, ['a', 'b', 'c'])'
    it will return:
        merge_all_subdicts(merge_all_subdicts(merge_all_subdicts(input_ensemble[a])[b])[c])

    the first element of the key_list will be the one that's merged first
    """

    key = key_list.pop(-1)
    if not key_list:
        # return "merge_all_subdicts(input_ensemble[{}])".format(key)
        return merge_all_subdicts(input_ensemble[key])
    else:
        # return "merge_all_subdicts({}[{}])".format(recursively_merge_all_subdicts(input_ensemble, key_list), key)
        return merge_all_subdicts(recursively_merge_all_subdicts(input_ensemble, key_list)[key])

def merge_all_subdicts(input_dict):
    """
    Function that merges all the subdicts of the given dict (not recursively, so first level only)
    """

    merged_results = {}
    for key in input_dict.keys():
        merged_results = merge_dicts(merged_results, input_dict[key])

    return merged_results

def merge_dicts(dict_1: dict, dict_2: dict) -> dict:
    """
    Function to merge two dictionaries in a custom way.

    Behaviour:

    When dict keys are only present in one of either:
        - we just add the content to the new dict

    When dict keys are present in both, we decide based on the value types how to combine them:
        - dictionaries will be merged by calling recursively calling this function again
        - numbers will be added
        - (opt) lists will be appended
        - In the case that the instances do not match: for now I will raise an error

    Args:
        dict_1: first dictionary
        dict_2: second dictionary

    Returns:
        Merged dictionary

    """

    # Set up new dict
    new_dict = collections.OrderedDict()  # TODO: check if this still necessary

    #
    keys_1 = dict_1.keys()
    keys_2 = dict_2.keys()

    # Find overlapping keys of both dicts
    overlapping_keys = set(keys_1).intersection(set(keys_2))

    # Find the keys that are unique
    unique_to_dict_1 = set(keys_1).difference(set(keys_2))
    unique_to_dict_2 = set(keys_2).difference(set(keys_1))

    # Add the unique keys to the new dict
    for key in unique_to_dict_1:
        # If these items are ints or floats, then just put them in
        if isinstance(dict_1[key], (float, int)):
            new_dict[key] = dict_1[key]
        # Else, to be safe we should deepcopy them
        else:
            copy_dict = dict_1[key]
            new_dict[key] = copy_dict

    for key in unique_to_dict_2:
        # If these items are ints or floats, then just put them in
        if isinstance(dict_2[key], (float, int)):
            new_dict[key] = dict_2[key]
        # Else, to be safe we should deepcopy them
        else:
            copy_dict = dict_2[key]
            new_dict[key] = copy_dict

    # Go over the common keys:
    for key in overlapping_keys:

        # If they keys are not the same, it depends on their type whether we still deal with them at all, or just raise an error
        if not type(dict_1[key]) is type(dict_2[key]):
            # Exceptions: numbers can be added
            if isinstance(dict_1[key], (int, float, np.float64)) and isinstance(
                dict_2[key], (int, float, np.float64)
            ):
                new_dict[key] = dict_1[key] + dict_2[key]

            # Exceptions: versions of dicts can be merged
            elif isinstance(
                dict_1[key], (dict, collections.OrderedDict)
            ) and isinstance(
                dict_2[key], (dict, collections.OrderedDict)
            ):
                new_dict[key] = merge_dicts(dict_1[key], dict_2[key])

            # If the above cases have not dealt with it, then we should raise an error
            else:
                print(
                    "Error key: {} value: {} type: {} and key: {} value: {} type: {} are not of the same type and cannot be merged".format(
                        key,
                        dict_1[key],
                        type(dict_1[key]),
                        key,
                        dict_2[key],
                        type(dict_2[key]),
                    )
                )
                raise ValueError

        # Here we check for the cases that we want to explicitly catch. Ints will be added,
        # floats will be added, lists will be appended (though that might change) and dicts will be
        # dealt with by calling this function again.
        else:
            # ints
            # Booleans (has to be the type Bool, not just a 0 or 1)
            if isinstance(dict_1[key], bool) and isinstance(dict_2[key], bool):
                new_dict[key] = dict_1[key] or dict_2[key]

            elif isinstance(dict_1[key], int) and isinstance(dict_2[key], int):
                new_dict[key] = dict_1[key] + dict_2[key]

            # floats
            elif isinstance(dict_1[key], float) and isinstance(dict_2[key], float):
                new_dict[key] = dict_1[key] + dict_2[key]

            # lists
            elif isinstance(dict_1[key], list) and isinstance(dict_2[key], list):
                new_dict[key] = dict_1[key] + dict_2[key]

            # Astropy quantities (using a dummy type representing the numpy array)
            elif isinstance(dict_1[key], type(np.array([1]) * u.m)) and isinstance(
                dict_2[key], type(np.array([1]) * u.m)
            ):
                new_dict[key] = dict_1[key] + dict_2[key]

            # dicts
            elif isinstance(dict_1[key], dict) and isinstance(dict_2[key], dict):
                new_dict[key] = merge_dicts(dict_1[key], dict_2[key])

            else:
                print(
                    "Object types {}: {} ({}), {} ({}) not supported.".format(
                        key,
                        dict_1[key],
                        type(dict_1[key]),
                        dict_2[key],
                        type(dict_2[key]),
                    )
                )
                raise ValueError

    #
    return new_dict

def merge_and_recurse_two_levels(result_dict, named_subkey_1, named_subkey_2):
    """
    Function to merge all subdicts, until it found the first level, then does the same until it found the second level. then gets recursive sum. The order of subkey_1 and subkey_2 must be the natural order they occur in in the ensemble

    TODO: Change this to a more general function where the named subkeys are a list.
    """

    # Make sure we loop until we find the first named subkey
    result_dict = merge_recursive_until_stop(result_dict, named_subkey_1)

    # Make sure that for each value key in the first named subkey dict we merge it until we found the correct named subkey 2:
    for value_subkey_1 in result_dict[named_subkey_1].keys():
        # replace the dict with the merged second subkey dict
        result_dict[named_subkey_1][value_subkey_1] = merge_recursive_until_stop(result_dict[named_subkey_1][value_subkey_1], named_subkey_2)

    # Then the levels below the second subkey should be recursive-summed
    result_dict = get_recursive_sum_for_subkeys_2d(result_dict, named_subkey_1, named_subkey_2)

    return result_dict

def get_recursive_sums_for_all_parameter(ensemble_dict, stopping_key='q_acc_don'):
    """
    Function to continuously find the recursive sum of the value keys of a parameter key. And then merge all subdicts to go a level deeper. 
    """

    #
    result_dict = OrderedDict()

    #
    calculating = True
    while calculating:
        current_parameter_dict = {}

        # Check if we only have 1 parameter key at this level
        if not len(list(ensemble_dict.keys())) == 1:
            print(ensemble_dict.keys())
            raise ValueError('Something went wrong: {}'.format(ensemble_dict.keys()))
        parameter_key = list(ensemble_dict.keys())[0]

        # Go over the value keys, calculate their total sum and add that to the new dict
        for subkey in ensemble_dict[parameter_key]:
            if not isinstance(ensemble_dict[parameter_key][subkey], (dict, OrderedDict)):
                current_parameter_dict[subkey] = ensemble_dict[parameter_key][subkey]
            else:
                current_parameter_dict[subkey] = get_recursive_sum(ensemble_dict[parameter_key][subkey])

        # add current dict to the total dict
        result_dict[parameter_key] = current_parameter_dict

        # Check if we are at the last one
        if parameter_key == stopping_key:
            calculating = False

        if calculating:
            ensemble_dict = merge_all_subdicts(ensemble_dict[parameter_key])

    #
    return result_dict


def find_columnames_recursively(ensemble_data, columnnames=None):
    """
    Function to find all the column names recursively

    This function should only be called on ensemble datasets that do not have differnt tree structures in them

    The first layer should be a namelayer
    """

    # get the column name
    if columnnames is None:
        new_columnnames = [list(ensemble_data.keys())[0]]
    else:
        new_columnnames = columnnames + [list(ensemble_data.keys())[0]]

    # Check if we are in the lowest layer
    next_layer_keys = list(ensemble_data[new_columnnames[-1]].keys())

    next_next_layer = ensemble_data[new_columnnames[-1]][next_layer_keys[0]]

    # Call itself or return if
    if isinstance(next_next_layer, (dict, OrderedDict)):
        return find_columnames_recursively(next_next_layer, columnnames=new_columnnames)

    return new_columnnames

def inflate_ensemble(ensemble_data):
    """
    Function to inflate an ensmeble, taking all the values for each datalayer and making a rectangular grid for it

    The first value should be a namelayer
    """

    parameter_name = list(ensemble_data.keys())[0]

    next_layer_keys = list(ensemble_data[parameter_name].keys())

    first_of_next_next_layer = ensemble_data[parameter_name][next_layer_keys[0]]
    is_final_layer = not isinstance(first_of_next_next_layer, (dict, OrderedDict))

    # if this is the final layer, then handle this layer with the dedicated function
    if is_final_layer:
        flattened_data = flatten_data_ensemble1d(ensemble_data, parameter_name).T
        return flattened_data

    # If its not the final layer, we should call this function again and return the result
    combined_array = None
    for valuekey in next_layer_keys:
        # Check if its the final layer
        next_next_layer = ensemble_data[parameter_name][valuekey]

        # combine result with the current value key
        res = inflate_ensemble(next_next_layer)

        # look at how many rows we have
        rows = res.shape[0]
        cols = res.shape[-1]

        try:
            converted_key = float(valuekey)
        except:
            converted_key = str(valuekey)

        # Get a column with the current valuekeys
        new_column = np.array([converted_key] * rows)

        # Create a new array with the shape that can fit the new column
        new_array = np.zeros((rows, cols + 1))

        # Set the new column as the first one and the rest in the rest
        new_array[:, 0] = new_column
        new_array[:, 1:] = res

        # Combine the arrays
        if combined_array is None:
            combined_array = new_array
        else:
            combined_array = np.append(combined_array, new_array, axis=0)

    # Return the results
    return combined_array


def inflate_ensemble_with_dataframes(ensemble_data, final_layer_name='yield_per_solarmass'):
    """
    Function to inflate an ensemble by using dataframes. taking all the values for each datalayer and making a rectangular grid for it

    This method allows string-based value-keys, whereas inflate_ensemble allows only value-keys that can be transformed into floats.

    The first value should be a namelayer
    """

    parameter_name = list(ensemble_data.keys())[0]

    next_layer_keys = list(ensemble_data[parameter_name].keys())

    first_of_next_next_layer = ensemble_data[parameter_name][next_layer_keys[0]]
    is_final_layer = not isinstance(first_of_next_next_layer, (dict, OrderedDict))

    # if this is the final layer, then handle this layer with the dedicated function
    if is_final_layer:
        flattened_data = flatten_data_ensemble1d(ensemble_data, parameter_name).T
        flattened_data_df = pd.DataFrame(flattened_data, columns=[parameter_name, final_layer_name])
        flattened_data_df[final_layer_name] = flattened_data_df[final_layer_name].astype(float)

        return flattened_data_df

    # If its not the final layer, we should call this function again and return the result
    combined_df = pd.DataFrame()
    for valuekey in next_layer_keys:
        # Check if its the final layer
        next_next_layer = ensemble_data[parameter_name][valuekey]

        # combine result with the current value key
        res = inflate_ensemble_with_dataframes(next_next_layer)
        cols_list = [parameter_name] + res.columns.tolist()

        res[parameter_name] = valuekey
        res = res[cols_list]

        # Combine dataframes
        combined_df = pd.concat([combined_df, res]).reset_index(drop=True)

    # Return the results
    return combined_df


def inflate_ensemble_with_lists(ensemble_data):
    """
    Function to inflate an ensemble by using dataframes. taking all the values for each datalayer and making a rectangular grid for it

    The first value should be a namelayer
    """
    parameter_name = list(ensemble_data.keys())[0]

    next_layer_keys = list(ensemble_data[parameter_name].keys())

    first_of_next_next_layer = ensemble_data[parameter_name][next_layer_keys[0]]
    is_final_layer = not isinstance(first_of_next_next_layer, (dict, OrderedDict))

    # if this is the final layer, then handle this layer with the dedicated function
    if is_final_layer:
        flattened_data = flatten_data_ensemble1d(ensemble_data, parameter_name)

        flattened_data_list = [[], []]
        flattened_data_list[0] = flattened_data[0].tolist()
        flattened_data_list[1] = flattened_data[1].astype(float).tolist()

        return flattened_data_list

    # If its not the final layer, we should call this function again and return the result
    combined_list = None
    for valuekey in next_layer_keys:
        next_next_layer = ensemble_data[parameter_name][valuekey]

        # combine result with the current value key
        res = inflate_ensemble_with_lists(next_next_layer)

        if combined_list is None:
            combined_list = [[] for _ in range(len(res) + 1)]
        valuekey_list = [valuekey for _ in range(len(res[0]))]

        # Create add the current key to the combined_list
        combined_list[0] = combined_list[0] + valuekey_list

        # add the results of the previous inflation to the results
        for res_list_i, res_list in enumerate(res):
            combined_list[res_list_i + 1] = combined_list[res_list_i + 1] + res_list

    # Return the results
    return combined_list













# reduced_ensemble_data = {
#     'a': {
#         '5': {
#             'b': {
#                 '1': 0.5,
#                 '2': 0.5,
#             }
#         },
#         '6': {
#             'b': {
#                 '3': 0.5,
#                 '4': 0.5,
#             }
#         },
#         '7': {
#             'b': {
#                 '8': 0.5,
#                 '9': 0.5,
#             }
#         }
#     }
# }

# reduced_ensemble_data = {
#     'z': {
#         '10': {
#             'a': {
#                 '5': {
#                     'b': {
#                         '1': 0.5,
#                         '2': 0.5,
#                     }
#                 },
#                 '6': {
#                     'b': {
#                         '3': 0.5,
#                         '4': 0.5,
#                     }
#                 },
#                 '7': {
#                     'b': {
#                         '8': 0.5,
#                         '9': 0.5,
#                         }
#                     }
#                 }
#             },
#         '11': {
#             'a': {
#                 '5': {
#                     'b': {
#                         '1': 0.5,
#                         '2': 0.5,
#                     }
#                 },
#                 '6': {
#                     'b': {
#                         '3': 0.5,
#                         '4': 0.5,
#                     }
#                 },
#                 '7': {
#                     'b': {
#                         '8': 0.5,
#                         '9': 0.5,
#                         }
#                     }
#                 }
#             },
#         }
#     }

