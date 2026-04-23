"""Packages containing power spectral analyses

"""

__all__ = ['perievent_histogram', 'tide']

def __getattr__(name):
    if name == 'perievent_histogram':
        from .peth import perievent_histogram
        return perievent_histogram
    if name == 'tide':
        from .tide import tide
        return tide
    raise AttributeError(f"module 'seapipe.stats' has no attribute '{name}'")
