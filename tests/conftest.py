import os
import pytest

TESTS_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(TESTS_DIR, 'data')


def data_path(filename):
    return os.path.join(DATA_DIR, filename)


def read_file(path):
    with open(path) as f:
        return f.read()
