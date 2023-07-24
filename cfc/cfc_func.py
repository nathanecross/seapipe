# -*- coding: utf-8 -*-
"""
Created on Mon May 13 13:10:09 2019

CFC Functions

@author: jordan
"""

from numpy import (any, arange, arctan2, argmin, around, array, asarray, 
                   atleast_2d, ceil, concatenate, cumsum, dot, empty, exp, floor, genfromtxt, 
                   hstack, hypot, identity, linalg, linspace, log, logical_and, 
                   max, mean, minimum, multiply, nan, nanmean, nanstd, ones, 
                   pi, power, prod, reshape, roll, sin, sqrt, squeeze, std, sum, 
                   transpose, vstack, where, zeros)
from os import listdir, mkdir, path
from scipy.stats import f as fdist
from scipy.special import lpn
import matplotlib.pyplot as plt
from pingouin import circ_r   
from safepickle import load
import warnings


def circ_wwtest(alpha1, alpha2, w1, w2, warnings=True):
    """Parametric Watson-Williams multi-sample test for equal means. Can be 
    used as a one-way ANOVA test for cicular data. Adapted for binned data."""        
    n1 = sum(w1)
    n2 = sum(w2)
    N = n1 + n2
    
    r1 = circ_r(alpha1, w1)    
    r2 = circ_r(alpha2, w2)
    
    r = circ_r(concatenate((alpha1, alpha2)), concatenate((w1, w2)))
    rw = sum((n1 * r1, n2 * r2)) / N
    
    if warnings:
        # check assumptions
        if N >= 11 and rw < .45:
            print('Warning: Test not applicable.'
                  'Average resultant vector length < 0.45.')
        elif N < 11 and N >= 7 and rw < .5:
            print('Test not applicable. Average number of samples per population '
                  '6 < x < 11 and average resultant vector length < 0.5.')
        elif N >= 5 and N < 7 and rw < .55:
            print('Test not applicable. Average number of samples per population '
                  '4 < x < 7 and average resultant vector length < 0.55.')
        elif N < 5:
            print('Test not applicable. '
                  'Average number of samples per population < 5.')
    
    # test statistic
    kk = circ_kappa(rw)
    beta = 1 + 3 / (8 * kk) # correction factor
    A = sum((n1 * r1, n2 * r2)) - r * N
    B = N - sum((n1 * r1, n2 * r2))
    
    F = beta * (N - 2) *  A / B
    pval = 1 - fdist.cdf(F, 1, N - 2)
    
    return F, pval
        
def circ_kappa(alpha, w=None):
    """ Computes an approximation to the ML estimate of the concentration 
    parameter kappa of the von Mises distribution."""    
    if isinstance(alpha, float):
        r = alpha
        N = 1        
    else:
        r = circ_r(alpha, w)
        N = len(alpha)
        
    if r < .53:
        kappa = (2 * r) + (r ** 3) + (5 * r ** (5/6))
    elif r >= .53 and r < .85:
        kappa = -.4 + 1.39 * r + .43 / (1 - r)
    else:
        kappa = 1 / (r ** 3 - 4 * r ** 2 + 3 * r)
        
    if N < 15 and N > 1:
        if kappa < 2:
            kappa = max(kappa - 2 / (N * kappa), 0)
        else:
            kappa = (N - 1) ** 3 * kappa / (N ** 3 + N)
            
    return kappa

def polar_plot(w, bin_edges):
    nbins = len(bin_edges) - 1
    width = 2 * pi / nbins
    f, ax = plt.subplots(1, 1, subplot_kw=dict(projection='polar'))
    ax.bar(bin_edges[:nbins], w, width=width, bottom=0.0)
    xL = ['0', '', '', '', r'+/- $\pi$', '', '', '']
    ax.set_xticklabels(xL)
    ax.set_yticks([])
    f.show()

def mean_amp(phase, amp, nbins=18):
    # Meanamp taking Hilbert transformed phase and amplitude signals (for comodulogram)
    # phase is the phase time series, type = array
    # amp is the amplitude time series, type = array
    # position is the array of phase bin left boundaries, type = array   
    # First, we break 2pi into phase bins
    # 0 degrees corresponds to the negative trough of the Phase freq
    # 180 degrees corresponds to the positive peak of the Phase freq       
    position = zeros(nbins)
    winsize = 2 * pi / nbins
    for i in range(nbins):
        position[i] = -pi + i * winsize
    position = roll(position, int(floor(nbins/4)-1), axis=-1)
    # Then, we compute mean amplitude in each phase bin:
    ma = zeros(nbins)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for i, j in enumerate(position):
            ma[i] = nanmean(amp[where(
                logical_and(phase < (j + winsize), phase >= j)
            )])        
    return ma

def klentropy(MeanAmp, axis=-1):    
    # Quantify the amount of amp modulation by means of a normalized 
    # Kullback-Leibler entropy index
    nbin = MeanAmp.shape[-1]
    MI = (log(nbin) - (-sum((MeanAmp / sum(MeanAmp, axis=axis, keepdims=True)) * \
        log((MeanAmp / sum(MeanAmp, axis=axis, keepdims=True))), axis=axis, 
        keepdims=True))) / log(nbin)
    return MI.squeeze()

# =============================================================================
# ALL NIGHT AMPBIN WITH STD DEV
# def _allnight_ampbin(a, nodat, nbins, norm=True, sd=False):
#     out = zeros((*a.shape[:-2], nbins))
#     out_sd = zeros((*a.shape[:-2], nbins))
#     for i in range(a.shape[0]): # part
#         for j in range(a.shape[1]): # night
#             allcyc = [x for x in a[i, j, 0, :] if x is not nodat]
#             allcyc = concatenate(allcyc)
#             if norm:
#                 allcyc /= allcyc.sum(axis=1, keepdims=True)
#             out[i, j, :] = nanmean(allcyc, axis=0)
#             if std:
#                 out_sd[i, j, :] = nanstd(allcyc, axis=0)
#                 out = (out, out_sd)                
#     return out
#
# OLD _allnight_ampbin
# def _allnight_ampbin(a, nodat, nbins, norm=True):
#     out = zeros((*a.shape[:-2], nbins))
#     for i in range(a.shape[0]): # part
#         for j in range(a.shape[1]): # night
#             allcyc = [x for x in a[i, j, 0, :] if x is not nodat]
#             allcyc = concatenate(allcyc)
#             if norm:
#                 allcyc /= allcyc.sum(axis=1, keepdims=True)
#             out[i, j, :] = nanmean(allcyc, axis=0)                
#     return out
# =============================================================================

def _allnight_ampbin(a, nodat, nbins, norm=True):
    out = zeros((*a.shape[:-2], nbins))
    for i in range(a.shape[0]): # part
            allcyc = [x for x in a[i] if x is not nodat]
            allcyc = concatenate(allcyc)
            if norm:
                allcyc /= allcyc.sum(axis=1, keepdims=True)
            out[i, :] = nanmean(allcyc, axis=0)                
    return out

def ab_loader(ab_file, nbins, shift, norm=True):
    with open(ab_file, 'rb') as f:
        ab = load(f)
    try:
        ab = _allnight_ampbin(ab, 0, nbins, norm=norm)
    except(ValueError):
        ab = ab / ab.sum(-1, keepdims=True)
        ab = ab.squeeze()    
    ab = roll(ab, shift, axis=-1)
    
    return ab

def _idx_to_phase(idx, vecbin):
    phase = zeros(idx.shape, dtype='O')
    for i in range(idx.shape[0]):
        for j in range(idx.shape[1]):
            one_idx = asarray(idx[i][j], dtype='int')
            phase[i][j] = vecbin[one_idx]
            
    return phase

def logsd(m, sd):
    print('M: ' + str(around(exp(m), 3)))
    print('SD: ' + str(around(sqrt((exp(sd ** 2) - 1) * exp(2 * m + sd ** 2)), 3)))
    
def logse(m, sd, N=15):
    print('M: ' + str(around(exp(m), 3)))
    print('SE: ' + str(around(
            sqrt((exp(sd ** 2) - 1) * exp(2 * m + sd ** 2)) / sqrt(N), 3)))

def cohend(x, y):
    mean_diff = x.mean() - y.mean()
    pooled_sd = sqrt((x.std() ** 2 + y.std() ** 2) / 2)
    return mean_diff / pooled_sd

    
def paplot(ma):
    f, axarr = plt.subplots(nrows=2, sharex=True)
    width = 2 * pi / len(ma)
    pos = arange(0, 4 * pi, width)
    x = linspace(0, 4 * pi, 1000)
    pa = axarr[0]
    pa.bar(pos, hstack((ma, ma)), width=width, color='k', edgecolor='w')
    pa.set_xticks([0, pi/2, pi, 3 * pi / 2, 2 * pi, 5 * pi / 2, 3 * pi, 7 * pi / 2, 4 * pi])
    pa.set_xticklabels(['0'] + [r'$\pi$/2', r'$\pi$', r'3$\pi$/2', r'2$\pi$'] * 2)
    sine = axarr[1]
    sine.plot(sin(x))
    sine.set_xticks([])
    sine.set_yticks([])
    plt.show()
    