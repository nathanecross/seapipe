"""
SEA main module
"""

__author__ = 'Nathan Cross'
__credits__ = 'Concordia University (Canada), University of Sydney (Australia)'


from os import path

here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'VERSION')) as f:
    __version__ = f.read().strip()


