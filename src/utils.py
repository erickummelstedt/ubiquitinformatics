import json

def match_assertion_error_contains(error_message, expected_parts):
    """
    Returns True if all strings in `expected_parts` are found in `error_message`.
    """
    return all(part in error_message for part in expected_parts)

def all_strings_exist(substrings, error_string):
    """
    Checks if all strings in the list `substrings` are present in `error_string`.

    :param substrings: List of strings to check for
    :param error_string: The string to search within
    :return: True if all substrings are found, False otherwise
    """
    return all(sub in error_string for sub in substrings)

def all_strings_exist_in_list(expected_substrings, error_strings):
    """
    Checks if all strings in each expected_substring exist in the corresponding string in error_strings.

    Args:
        error_strings (List[str]): A list of error_string strings.
        expected_substrings (List[List[str]]): A list of expected_substrings of substrings to check against each error_string string.

    Returns:
        List[bool]: A list of booleans indicating match success for each error_string/expected_substring pair.
    """
    if len(error_strings) != len(expected_substrings):
        raise ValueError("Length of error_string and expected_substrings must be equal.")

    result = []
    for error_string, expected_substring in zip(error_strings, expected_substrings):
        match = all(sub in error_string for sub in expected_substring)
        result.append(match)

    return result

def convert_json_to_dict(parent_dictionary):
    """
    Converts a JSON string to a dictionary if necessary.
    If the input is already a dictionary, it remains unchanged.
    Raises a ValueError if the JSON format is invalid.

    :param input_data: JSON string or dictionary
    :return: Dictionary representation of the input data
    """

    if isinstance(parent_dictionary, str):
        try:
            # Ensure correct JSON format by replacing single quotes with double quotes
            formatted_json = parent_dictionary.replace("'", "\"")
            return json.loads(formatted_json)
        except json.JSONDecodeError as e:
            raise ValueError("Invalid JSON format: Unable to parse the string") from e
    elif isinstance(parent_dictionary, dict):
        return parent_dictionary
    else:
        raise TypeError("Input must be a dictionary or a JSON string")
