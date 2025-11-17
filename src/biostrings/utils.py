import biocutils as ut

__author__ = "Jayaram Kancherla"
__copyright__ = "Jayaram Kancherla"
__license__ = "MIT"


def _sanitize_metadata(metadata):
    if metadata is None:
        return {}
    elif not isinstance(metadata, dict):
        metadata = dict(metadata)

    return metadata


def _sanitize_names(names):
    if names is None:
        return None
    elif not isinstance(names, list):
        names = ut.Names(names)

    return names
