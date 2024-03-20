"""Packages containing power spectral analyses

"""
from .psa import (Spectrum, default_epoch_opts, default_event_opts, 
                  default_fooof_opts, default_filter_opts, default_frequency_opts, 
                  default_general_opts)
from .spectrogram import event_spectrogram, event_spectrogram_grouplevel

