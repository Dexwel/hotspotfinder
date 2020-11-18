import gzip


def file_open(file, mode='wt'):
    """Open file. Handles normal and gzip files"""
    if file.endswith('.gz'):
        open_function = gzip.open
    else:
        open_function = open
    return open(file, mode=mode)
