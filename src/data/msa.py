"""Module which contains helper functions for storing/retrieving data related to multi sequence
alignments generated using hmmalign
"""

def get_aligned_bgc_list(database):
    """Returns a list of ids of genomes which already have predicted bgcs"""
    predicted_bgcs = [
        row["bgc_id"] for row in database.select(
            "bgc_status",
            "where status > 3",
            props=["bgc_id"])]
    return predicted_bgcs
