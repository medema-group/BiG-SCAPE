"""Contains code to calculate the DSS for a pair of BGCs. Also contains a helper
function to perform calculation for all BGCs in a bin using just-in-time execution of
HHMAlign
"""

# from python
# from dependencies
# from other modules

# from this module


def get_aligned_string_dist(string_a: str, string_b: str) -> float:
    """Calculate a simple distance between two strings of equal length from an MSA

    Strings must be equal lengths.

    Args:
        string_a (str): String to calculate distance for
        string_b (str): String to calculate distance for

    Raises:
        ValueError: Raised when string lengths do not match

    Returns:
        float: Simple distance of the two passed strings
    """
    if len(string_a) != len(string_b):
        raise ValueError(
            "String A and String B length difference in get_aligned_string_dist"
        )

    gaps = 0
    matches = 0

    for char_idx in range(len(string_a)):
        if string_a[char_idx] == string_b[char_idx]:
            if string_a[char_idx] == "-":
                gaps += 1
            else:
                matches += 1

    similarity = matches / (len(string_a) - gaps)

    return 1 - similarity
