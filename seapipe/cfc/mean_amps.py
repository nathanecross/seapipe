# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 15:51:24 2021

@author: Nathan Cross
"""
from termcolor import colored
from os import listdir, mkdir, path, walk
from . cfc_func import _allnight_ampbin, circ_wwtest, mean_amp, klentropy
from utils.misc import bandpass_mne, laplacian_mne, notch_mne, notch_mne2
import copy
from math import degrees, radians
import mne
import matplotlib.pyplot as plt
from numpy import (angle, append, argmax, array, arange, asarray, ceil, concatenate, 
                   empty, histogram, interp, isnan, linspace, log, logical_and, mean, 
                   median, nan, nanmean, ndarray, newaxis, ones, pi, random, repeat, 
                   reshape, roll, sin, size, squeeze, sqrt, std, sum, tile, where, zeros) 
from numpy.matlib import repmat
from pandas import DataFrame, concat, read_csv
from pathlib import Path
from safepickle import dump, load
from pingouin import (circ_mean, circ_r, circ_rayleigh, circ_corrcc, circ_corrcl)
from scipy.signal import hilbert
from scipy.stats import zscore
import sys
from tensorpac import Pac, EventRelatedPac
from wonambi import Dataset
from wonambi.trans import fetch
from wonambi.attr import Annotations 
from wonambi.detect.spindle import transform_signal


def pac_it(rec_dir, xml_dir, out_dir, part, visit, cycle_idx, chan, rater, stage,
               polar, grp_name, cat, evt_type, buffer, ref_chan, nbins, idpac, 
               fpha, famp, dcomplex, filtcycle, width, min_dur, band_pairs,
               adap_bands=(False,False),
               filter_opts={'notch':False,'notch_harmonics':False, 'notch_freq':None,
                            'laplacian':False, 'lapchan':None,'laplacian_rename':False, 
                            'oREF':None,'chan_rename':False,'renames':None},
               progress=True):

    '''
    This script runs Phase Amplitude Coupling analyses on sleep EEG data. 
    The method for calculating PAC is set by the parameter <idpac>. 
    For more information on the available methods, refer to the documentation of 
    tensorpac (https://etiennecmb.github.io/tensorpac/index.html) or the article
    (Combrisson et al. 2020, PLoS Comp Bio: https://doi.org/10.1371/journal.pcbi.1008302)
    
    The script does the following:
        1. Extracts the EEG signal at each event specified by <evt_type> ± a buffer
           on either side of length (in sec) specified by <buffer>.
        2. For these EEG segments, filters the signal within a given frequency range
           specified by <fpha> to obtain the phase, and again within a given frequency 
           range specified by <famp>.  
        3. FOR EACH EACH EVENT: the instantaneous amplitude of the signal filtered 
           within a given frequency range specified by <famp> will be calculated 
           via the Hilbert transform, and the amplitude will be averaged across a 
           set number of phase bins specified by <nbins>. The phase bin with the 
           maxmimum mean amplitude will be stored.
        4. ACROSS ALL EVENTS: the average phase bin with the maximum amplitude will 
           be calculated (circular mean direction).
        5. The filtered events will also be concatenated and stored in blocks of 50,
           so that the PAC strength (method pecficied by the 1st entry in <idpac>)
           can be calculated AND surrogates can be accurately generated to test for
           the significance of PAC in each participant. The generation of surrogates
           is specified by the 2nd entry in <idpac>, and the correction of PAC 
           strength is also calculated, specified by the 3rd entry in <idpac>.
        6. Other metrics are also calculated for each participant and visit, notably:
            - mean vector length (given from mean circular calculation)
            - correlation between amplitudes (averaged over all events) and the phase 
               giving sine wave
            - Rayleigh test for non-uniformity of circular data (sig. test for 
                                                                 preferred phase)
            
           
    If laplacian = True then a Laplacian spatial filter will be applied to remove high frequency EMG 
    noise. In this scenario you will need to provide a list of channel names to include in the laplacian
    spatial filtering. 
                    ## WARNING: TEST WHAT THE LAPLACIAN FILTER DOES TO THE POWER SPECTRUM BEFORE USING
                                THIS OPTION
    
    If adap_bands = (True,True) then the (phase,amplitude) signal will be filtered within an adapted 
    frequency range for each individual subject or recording.
    
    The output provided by this script will be an array of size:
        [#cycles x #bins] (if cycle_idx is not None)
        [1 x #nbins]      (if cycle_idx is None)
                                        - corresponding to the mean amplitude of the signal 
                                         (across all cycles or events) per phase bin.

    '''
    
    
    
    
    # Make output directory
    if not path.exists(out_dir):
        mkdir(out_dir)
    
    
    ## BIDS CHECKING
    # Check input participants
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(rec_dir)
            part = [ p for p in part if not '.' in p]
    else:
        print('')
        print(colored('ERROR |', 'red', attrs=['bold']),
              colored("'part' must either be an array of subject ids or = 'all'" ,'cyan', attrs=['bold']))
        print('')
    
    # Check input visits
    if isinstance(visit, list):
        None
    elif visit == 'all':
        lenvis = set([len(next(walk(rec_dir + x))[1]) for x in part])
        if len(lenvis) > 1:
            print(colored('WARNING |', 'yellow', attrs=['bold']),
                  colored('number of visits are not the same for all subjects.',
                          'white', attrs=['bold']))
            print('')
        visit = list(set([y for x in part for y in listdir(rec_dir + x)  if '.' not in y]))
    else:
        print('')
        print(colored('ERROR |', 'red', attrs=['bold']),
              colored("'visit' must either be an array of subject ids or = 'visit'" ,
                      'cyan', attrs=['bold']))
        print('')
    
    # Loop through participants and visits
    
    part.sort()
    for i, p in enumerate(part):
        
        ###########################            DEBUGGING              ###########################
        with open(out_dir + f'/debug_{p}.txt', 'w') as f:
            f.write('making participant output directory')
        ###########################            DEBUGGING              ###########################
        
        if not path.exists(out_dir + '/' + p):
            mkdir(out_dir + '/' + p)
        
        if visit == 'all':
            visit = listdir(xml_dir + '/' + p)
            visit = [x for x in visit if not '.' in x]
        visit.sort()    
        for j, vis in enumerate(visit): 
            if not path.exists(xml_dir + '/' + p + '/' + vis + '/'):
                print(colored('WARNING |', 'yellow', attrs=['bold']),
                      colored(f'input folder missing for Subject {p}, visit {j} ({vis}), skipping...',
                              'white', attrs=['bold']))
                print('')
                continue
            else:
                
                ###########################            DEBUGGING              ###########################
                with open(out_dir + f'/debug_{p}.txt', 'w') as f:
                    f.write(f'making visit {j} output directory')
                ###########################            DEBUGGING              ###########################
                
                if not path.exists(out_dir + '/' + p + '/' + vis):
                    mkdir(out_dir + '/' + p + '/' + vis)
                rec_file = [s for s in listdir(rec_dir + '/' + p + '/' + vis) if 
                            (".edf") in s or ('.rec') in s or ('.eeg')  in s if not s.startswith(".")]
                xml_file = [x for x in listdir(xml_dir + '/' + p + '/' + vis) if 
                            x.endswith('.xml') if not x.startswith(".")] 
                
                # Open recording and annotations files (if existing)
                if len(xml_file) == 0:
                    print(colored('WARNING |', 'yellow', attrs=['bold']),
                          colored(f'annotations does not exist for Subject {p}, visit {j} ({vis}) - check this. Skipping...',
                                  'white', attrs=['bold']))
                    print('')
                elif len(xml_file) >1:
                    print(colored('WARNING |', 'yellow', attrs=['bold']),
                          colored(f'multiple annotations files exist for Subject {p}, visit {j} ({vis}) - check this. Skipping...',
                                  'white', attrs=['bold']))
                    print('')
                else:
                    
                    ###########################            DEBUGGING              ###########################
                    with open(out_dir + f'/debug_{p}.txt', 'w') as f:
                        f.write('opening participant edf and xml')
                    ###########################            DEBUGGING              ###########################
                    
                    dset = Dataset(rec_dir + '/' + p + '/' + vis + '/' + rec_file[0]) 
                    annot = Annotations(xml_dir + '/' + p + '/' + vis + '/' + xml_file[0], 
                                        rater_name=rater)
                    
                    # Get sleep cycles
                    if cycle_idx is not None and cat[0] == 1:
                        all_cycles = annot.get_cycles()
                        scycle = [all_cycles[i - 1] for i in cycle_idx if i <= len(all_cycles)]
                        
                        # Create output array (length #cycles)
                        all_ampbin = zeros((6), dtype='object')
                    else:
                         
                        # Create output array (length 1)
                        all_ampbin = zeros((1), dtype='object')
                        scycle = [None]
                    
                    # Loop through channels
                    for k, ch in enumerate(chan):
                        
                        print(f'Reading data for {p}, visit {vis}, channel {ch}')
                        chan_full = ch + ' (' + grp_name + ')'
                        
                        # Check for adapted bands
                        if adap_bands[0] is True:
                            f_pha = fpha[ch][p + '_' + vis]
                            print(f'Using adapted phase frequency for {p}: {round(f_pha[0],2)}-{round(f_pha[1],2)} Hz')
                        else:
                            f_pha = fpha
                        if adap_bands[1] is True:
                            f_amp = famp[ch][p + '_' + vis]
                            print(f'Using adapted amplitude frequency for {p}: {round(f_amp[0],2)}-{round(f_amp[1],2)} Hz')
                        else:
                            f_amp = famp
                            
                        # Loop through sleep cycles 
                        for l, cyc in enumerate(scycle):
                            print('')
                            # Select and read data
                            if cycle_idx is not None:
                                print(f'Analysing, cycle {l+1}')
                            else:
                                print('Analysing, all cycles')
                            print('')
                            print('Using filter settings:')
                            print('')
                            print(colored('Notch filter:','white', attrs=['bold']),
                                  colored(f"{filter_opts['notch']}", 'yellow', attrs=['bold']))
                            print(colored('Notch harmonics filter:','white', attrs=['bold']),
                                  colored(f"{filter_opts['notch_harmonics']}", 'yellow', attrs=['bold']))
                            print(colored('Laplacian filter:','white', attrs=['bold']),
                                  colored(f"{filter_opts['laplacian']}", 'yellow', attrs=['bold']))
                            
                            cat = list(cat)
                            cat[2] = 0 #enforce non-concatenation of signal
                            cat[3] = 0 #enforce non-concatenation of event types
                            cat = tuple(cat)
                            
                            segments = fetch(dset, annot, cat=cat, chan_full=[chan_full], 
                                             cycle=[cyc], evt_type=evt_type, stage=stage,
                                             buffer=buffer)
                            if filter_opts['laplacian'] or filter_opts['notch'] or filter_opts['notch_harmonics']:
                                chans = filter_opts['lapchan']
                            else:
                                chans = [ch]
                            segments.read_data(chan=chans, ref_chan=ref_chan)
                            
                            if len(segments) <1:
                                print(colored('WARNING |', 'yellow', attrs=['bold']),
                                      colored('No segments found.',
                                              'white', attrs=['bold']))
                                
                            
                            # Define output & PAC object
                            
                            pac = Pac(idpac=idpac, f_pha=f_pha, f_amp=f_amp, dcomplex=dcomplex, 
                                      cycle=filtcycle, width=width, n_bins=nbins)
                            
                            
                            ###########################            DEBUGGING              ###########################
                            with open(out_dir + f'/debug_{p}.txt', 'w') as f:
                                f.write('Computing PAC')
                            ###########################            DEBUGGING              ###########################
                            
                            nsegs=[]
                            
                            if cat[1] == 0:
                                print('Splitting stages')
                                for s, st in enumerate(stage):
                                    segs = [s for s in segments if st in s['stage']]
                                    nsegs.append(segs)
                                    
                            else:
                                nsegs = [segments]
                            
                            for s in range(len(nsegs)):
                                segments = nsegs[s]
                                ampbin = zeros((len(segments), nbins))
                                
                                if cat[1] == 1:
                                    stagename = ''.join(stage) 
                                else: 
                                    stagename = stage[s]
                                    
                                print('')
                                print(f'Stage {stagename}')
                                print('Calculating mean amplitudes')
                                print(f'No. Segments = {len(segments)}')
                               
                                # Create blocks
                                ms = int(ceil(len(segments)/50))
                                longamp = zeros((ms,50),dtype=object) # initialise (blocked) ampltidue series
                                longpha = zeros((ms,50),dtype=object) # initialise (blocked) phase series 
                                
                                z=0
                                for m, seg in enumerate(segments):
                                    
                                    # Print out progress
                                    if progress:
                                        z +=1
                                        j = z/len(segments)
                                        sys.stdout.write('\r')
                                        sys.stdout.write(f"Progress: [{'=' * int(50 * j):{50}s}] {int(100 * j)}%")
                                        sys.stdout.flush()
                                    
                                    # Select data from segment
                                    data = seg['data']
                                    
                                    # Find sampling frequency
                                    s_freq = data.s_freq

                                    # Apply filtering (if necessary)
                                    if filter_opts['notch']:
                                        selectchans = list(data.chan[0])
                                        data.data[0] = notch_mne(data, oREF=filter_opts['oREF'], 
                                                                    channel=selectchans, 
                                                                    freq=filter_opts['notch_freq'],
                                                                    rename=filter_opts['chan_rename'],
                                                                    renames=filter_opts['renames'])
                                    
                                    if filter_opts['notch_harmonics']: 
                                        selectchans = list(data.chan[0])
                                        data.data[0] = notch_mne2(data, oREF=filter_opts['oREF'], 
                                                                  channel=selectchans,
                                                                  rename=filter_opts['chan_rename'],
                                                                  renames=filter_opts['renames'])
                                    
                                    if filter_opts['laplacian']:
                                        data = laplacian_mne(data, oREF=filter_opts['oREF'], channel=ch, 
                                                             ref_chan=ref_chan, 
                                                             laplacian_rename=filter_opts['laplacian_rename'], 
                                                             renames=filter_opts['renames'])
                                        dat = data[0]
                                    else:
                                        dat = data()[0][0]
                                    
                                    
                                    # Check polarity of recording
                                    if isinstance(polar, list):
                                        polarity = polar[i]
                                    else:
                                        polarity = polar
                                    if polarity == 'opposite':
                                        dat = dat*-1 
                                        
                                    # Obtain phase signal
                                    pha = squeeze(pac.filter(s_freq, dat, ftype='phase'))
                                    
                                    if len(pha.shape)>2:
                                        pha = squeeze(pha)
                                    
                                    # Obtain amplitude signal
                                    amp = squeeze(pac.filter(s_freq, dat, ftype='amplitude'))
                                    if len(amp.shape)>2:
                                        amp = squeeze(amp)
                                    
                                    # Extract signal (minus buffer)
                                    nbuff = int(buffer * s_freq)
                                    minlen = s_freq * min_dur
                                    if len(pha) >= 2 * nbuff + minlen:
                                        pha = pha[nbuff:-nbuff]
                                        amp = amp[nbuff:-nbuff]
                                                                   
                                        
                                    # Apply phase correction for hilbert transform
                                    #pha = roll(pha, int(pi/2*s_freq), axis=-1)
                                    
                                    
                                    # Put data in blocks (for surrogate testing)
                                    longpha[m//50, m%50] = pha
                                    longamp[m//50, m%50] = amp
                                    
                                    # Put data in long format (for preferred phase)
                                    ampbin[m, :] = mean_amp(pha, amp, nbins=nbins)
                                    
                                
                                # if number of events not divisible by block length
                                # pad incomplete final block with randomly resampled events
                                rem = len(segments) % 50
                                if rem > 0:
                                    pads = 50 - rem
                                    for pad in range(pads):
                                        ran = random.randint(0,rem)
                                        longpha[-1,rem+pad] = longpha[-1,ran]
                                        longamp[-1,rem+pad] = longamp[-1,ran]
                                    
                                # Calculate Coupling Strength
                                methods = {1: 'Mean Vector Length (MVL) :cite:`canolty2006high`',
                                           2 : 'Modulation Index (MI) :cite:`tort2010measuring`',
                                           3 : 'Heights Ratio (HR) :cite:`lakatos2005oscillatory`',
                                           4 : 'ndPAC :cite:`ozkurt2012statistically`',
                                           5 : 'Phase-Locking Value (PLV) :cite:`penny2008testing,lachaux1999measuring`',
                                           6 : 'Gaussian Copula PAC (GCPAC) :cite:`ince2017statistical`'}
                                surrogates = {0 :' No surrogates', 
                                              1 : 'Swap phase / amplitude across trials :cite:`tort2010measuring`',
                                              2 : 'Swap amplitude time blocks :cite:`bahramisharif2013propagating`',
                                              3 : 'Time lag :cite:`canolty2006high`'}
                                corrections = {0 : 'No normalization',
                                               1 : 'Substract the mean of surrogates',
                                               2 : 'Divide by the mean of surrogates',
                                               3 : 'Substract then divide by the mean of surrogates',
                                               4 : 'Z-score'}
                                print('')
                                print('Calculating coupling strength.')
                                print(f'Using method {methods[idpac[0]]}.')
                                print(f'Surrogate method: {surrogates[idpac[1]]}.')
                                print('Correcting strength using method:')
                                print(f'{corrections[idpac[2]]}.')
                                print('')
                                mi = zeros((longamp.shape[0],1))
                                mi_pv = zeros((longamp.shape[0],1))
                                for row in range(longamp.shape[0]): 
                                    amp = zeros((1))   
                                    pha = zeros((1)) 
                                    for col in range(longamp.shape[1]):
                                        pha = concatenate((pha,longpha[row,col]))
                                        amp = concatenate((amp,longamp[row,col]))
                                    pha = reshape(pha,(1,1,len(pha)))
                                    amp = reshape(amp,(1,1,len(amp)))
                                    mi[row] = pac.fit(pha, amp, n_perm=400,random_state=5,
                                                 verbose=False)[0][0]
                                    mi_pv[row] = pac.infer_pvalues(p=0.95, mcp='fdr')[0][0]

                                ## Calculate preferred phase
                                print('Caclulating preferred phase.')
                                print('')
                                ampbin = ampbin / ampbin.sum(-1, keepdims=True) # normalise amplitude
                                ampbin = ampbin.squeeze()
                                ampbin = ampbin[~isnan(ampbin[:,0]),:] # remove nan trials
                                ab = ampbin
                                
                                # Create bins for preferred phase
                                vecbin = zeros(nbins)
                                width = 2 * pi / nbins
                                for n in range(nbins):
                                    vecbin[n] = n * width + width / 2  
                                
                                # Calculate mean direction (theta) & mean vector length (rad)
                                ab_pk = argmax(ab,axis=1)
                                theta = circ_mean(vecbin,histogram(ab_pk,bins=nbins, 
                                                                   range=(0,nbins))[0])
                                theta_deg = degrees(theta)
                                if theta_deg < 0:
                                    theta_deg += 360
                                rad = circ_r(vecbin, histogram(ab_pk,bins=nbins)[0], d=width)
                                
                                # Take mean across all segments/events
                                ma = nanmean(ab, axis=0)
                                
                                # Correlation between mean amplitudes and phase-giving sine wave
                                sine = sin(linspace(-pi, pi, nbins))
                                sine = interp(sine, (sine.min(), sine.max()), (ma.min(), ma.max()))
                                rho, pv1 = circ_corrcc(ma, sine)
        
                                # Rayleigh test for non-uniformity of circular data
                                ppha = vecbin[ab.argmax(axis=-1)]
                                z, pv2 = circ_rayleigh(ppha)
                                pv2 = round(pv2,5)
                                
                                # Prepare filename
                                if cat[0] == 1:
                                    cyclename = 'wholenight'
                                else: 
                                    cyclename = f'cycle{l+1}'


                                # Save cfc metrics to dataframe
                                d = DataFrame([mean(pac.pac), mean(mi), median(mi_pv), theta, 
                                               theta_deg, rad, rho, pv1, z, pv2])
                                d = d.transpose()
                                d.columns = ['mi','mi_norm','sig','pp_rad','ppdegrees','mvl',
                                             'rho','pval','rayl','pval2']
                                d.to_csv(path_or_buf=out_dir + '/' + p + '/' + vis + '/' + 
                                         p + '_' + vis + '_' + ch + '_' + stagename + '_' + 
                                          '_' + cyclename + '_' + band_pairs + 
                                          '_cfc_params.csv', sep=',')
                                
                                # Save binned amplitudes to pickle file
                                with open(out_dir + '/' + p + '/' + vis + '/' + p + '_' + 
                                          vis + '_' + ch + '_' + stagename + '_' + 
                                          '_' + cyclename + '_' + band_pairs + 
                                          '_mean_amps.p', 'wb') as f:
                                     dump(ab, f)
   
    print('The function pac_it completed without error.')                             


def pac_it_2(rec_dir, xml_dir, out_dir, part, visit, cycle_idx, chan, rater, stage,
               polar, grp_name, cat, target, probe, buffer, ref_chan, nbins, idpac, 
               fpha, famp, dcomplex, filtcycle, width, min_dur, band_pairs,
               adap_bands=(False,False),
               filter_opts={'notch':False,'notch_harmonics':False, 'notch_freq':None,
                            'laplacian':False, 'lapchan':None,'laplacian_rename':False, 
                            'oREF':None,'chan_rename':False,'renames':None},
               progress=True):
   
    '''
    This script runs Phase Amplitude Coupling analyses on sleep EEG data. 
    The function is essentially the same as the function pac_it (above), but with
    a notable exception that it runs exclusively on coupled events that co-occur in time,
    as specified by the argument <target> (the function event_synchrony needs to be run 
    first, see cfc.synchrony). 
    The function will extract these coupled (co-occuring) events (e.g. SO+spindle), and 
    calculate the amplitude of the event specified by the argument <probe> (which should 
    be one of the co-occuring events, e.g. spindle). The timing (location) of the peak 
    amplitude of the probe event is then calcualted relative to the start time of the 
    <target> event, and this difference is calculated in terms of relative phase bins from
    the start of the <target> event. Any amplitude peaks <0 (ie. before start of the event) 
    or >360 (ie. after the end of the event) - are discarded.
    The modulation index is also calculated on the coupled events specified by <target>, 
    for only those events where the amplitude peak <probe> is betweem 0-360 degrees.
    '''
    
    # Start function
    if not path.exists(f'{out_dir}/'):
        mkdir(f'{out_dir}/')
    
    ## BIDS CHECKING
    
    # Check input participants
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(rec_dir)
            part = [p for p in part if not '.' in p]
            part.sort()
    else:
        print('')
        print(colored('ERROR |', 'red', attrs=['bold']),
              colored("'part' must either be an array of subject ids or = 'all'" ,'cyan', attrs=['bold']))
        print('')
    
    # Check input visits
    if isinstance(visit, list):
        None
    elif visit == 'all':
        lenvis = set([len(next(walk(rec_dir + x))[1]) for x in part])
        if len(lenvis) > 1:
            print(colored('WARNING |', 'yellow', attrs=['bold']),
                  colored('number of visits are not the same for all subjects.',
                          'white', attrs=['bold']))
            print('')
        visit = list(set([y for x in part for y in listdir(rec_dir + x)  if '.' not in y]))
        visit.sort()
    else:
        print('')
        print(colored('ERROR |', 'red', attrs=['bold']),
              colored("'visit' must either be an array of subject ids or = 'visit'" ,
                      'cyan', attrs=['bold']))
        print('')
    
    # Run PAC
    for i, p in enumerate(part):
        if not path.exists(f'{out_dir}/{p}'):
            mkdir(f'{out_dir}/{p}')
        for v, vis in enumerate(visit):
            print('')
            print('')
            print(f'Opening participant {p}, visit {vis}')
            if not path.exists(f'{out_dir}/{p}/{vis}/'):
                mkdir(f'{out_dir}/{p}/{vis}/')
            out = zeros((1,10)) # intialise output dataframe
            for c, ch in enumerate(chan):
                print('')
                print(f'Analysing PAC for channel {ch}')
                chan_full = f'{ch} ({grp_name})'
                
                # Open edf and get information
                rec_file = [s for s in listdir(rec_dir + '/' + p + '/' + vis) if 
                            (".edf") in s or ('.rec') in s or ('.eeg')  in s if not s.startswith(".")]
                dataset = Dataset(f'{rec_dir}/{p}/{vis}/{rec_file[0]}')
                s_freq = dataset.header['s_freq']
                nbuff = int(buffer*s_freq)
                
                # Open xml file and get information
                xml_file = [x for x in listdir(f'{xml_dir}/{p}/{vis}/') if '.xml' in x]
                annot = Annotations(f'{xml_dir}/{p}/{vis}/{xml_file[0]}')
                
                
                # EXTRACT BOTH EVENT TYPES
                # Events = (lower freq) phase signal
                events = [x for x in annot.get_events(target) if chan_full in x['chan']]
                events = events[1:]
                evt_mask = ones((len(events)), dtype=bool)

                
                # Probes = (higher freq) phase signal
                probes = [x for x in annot.get_events(probe) if chan_full in x['chan']]
                segs = fetch(dataset, annot, cat=(0,0,0,0), evt_type=[probe], cycle=None, 
                                 chan_full=[chan_full], reject_epoch=True, buffer=buffer,
                                 reject_artf = ['Artefact', 'Arou', 'Arousal'])
                segs.read_data([ch], ref_chan, grp_name=grp_name)
                
                
                # FIND SYNCRONOUS EVENTS
                # Vectorize start and end times and set up for broadcasting
                ev_beg = asarray([x['start'] for x in events])[:, newaxis]
                ev_end = asarray([x['end'] for x in events])[:, newaxis]
                pb_beg = asarray([x['start'] for x in probes])[newaxis, :]
                pb_end = asarray([x['end'] for x in probes])[newaxis, :]
                

                # Subtract every end by every start and find overlaps
                ev_minus_pb = ev_end - pb_beg # array of shape (len(det), len(std))
                pb_minus_ev = pb_end - ev_beg    
                overlapping = logical_and(ev_minus_pb > 0, pb_minus_ev > 0)
                
                # Check for non-overlapping events and remove
                if sum(sum(overlapping,axis=1)) > 0: # failsafe if all events are not coupled
                    emptyrows = sum(overlapping,axis=1)!=0
                try:
                    idx = [where(overlapping[x,:]==True)[0][0] if emptyrows[x] else 0 
                           for x in range(0,len(overlapping)) ]
                except:
                    print('')
                    print(colored('ERROR |', 'red', attrs=['bold']),
                          colored(f"There is an event '{target}' with no synced event '{probe}'." ,
                                  'cyan', attrs=['bold']))
                    print(colored("Check event types 'target' & 'probe', or try rerunning event synchrony." ,
                                  'cyan', attrs=['bold']))
                    print('')
                    return
                
                ## CALCULATE PREFERRED PHASE
                # Define frequency bands for PAC
                if adap_bands[0] is True:
                    f_pha = fpha[ch][p + '_' + vis]
                    print(f'Using adapted bands for phase: {round(f_pha[0],2)}-{round(f_pha[1],2)} Hz')
                elif adap_bands[0] is False:
                    f_pha = fpha
                    print(f'Using fixed bands for phase: {round(f_pha[0],2)}-{round(f_pha[1],2)}')
                    
                if adap_bands[1] is True:
                    f_amp = famp[ch][p + '_' + vis]
                    print(f'Using adapted bands for ampltidue: {round(f_amp[0],2)}-{round(f_amp[1],2)} Hz')
                elif adap_bands[1] is False:
                    f_amp = famp
                    print(f'Using fixed bands for ampltidue: {round(f_amp[0],2)}-{round(f_amp[1],2)}')
                    
                print('')
                print('Extracting preferred phase of coupled events:')
                print(f'{target} <-> {probe}')
                
                print('')
                print('Using filter settings:')
                print('')
                print(colored('Notch filter:','white', attrs=['bold']),
                      colored(f"{filter_opts['notch']}", 'yellow', attrs=['bold']))
                print(colored('Notch harmonics filter:','white', attrs=['bold']),
                      colored(f"{filter_opts['notch_harmonics']}", 'yellow', attrs=['bold']))
                print(colored('Laplacian filter:','white', attrs=['bold']),
                      colored(f"{filter_opts['laplacian']}", 'yellow', attrs=['bold']))
                print('')
                
                # Define PAC object
                pac = Pac(idpac=idpac, f_pha=fpha, f_amp=f_amp, dcomplex='hilbert', 
                          cycle=(3,6), width=7, n_bins=nbins)
                
                # Prepare outputs
                pkbin = zeros((len(events),1))
                z=0 #(progress param)
                for e, ev in enumerate(events):
                    
                    # Print out progress
                    if progress:
                        z +=1
                        j = z/len(events)
                        sys.stdout.write('\r')
                        sys.stdout.write(f"Progress: [{'=' * int(50 * j):{50}s}] {int(100 * j)}%")
                        sys.stdout.flush()            
                    
                    # Get times of base events (e.g. SOs)
                    (evbeg,evend) = (ev['start'],ev['end'])
                    
                    # Get times of probes (e.g. spindles)
                    pb = probes[idx[e]]
                    (pbbeg,pbend) = (pb['start'],pb['end'])
    
                    # Extract amplitude signal
                    seg = segs[idx[e]]
                    data = seg['data']()[0][0]
                    amp = squeeze(pac.filter(s_freq, data, ftype='amplitude'))
                    amp = amp[nbuff:-nbuff]
                    
                    # Find time of peak
                    peak = argmax(amp)                          # index of pk amplitude in event
                    pk_abs = peak/s_freq + pbbeg                # time of pk relative to start of recording
                    pk_rel = (pk_abs - evbeg)*s_freq            # time of pk relative to start of base event
                    window = int((evend - evbeg)*s_freq)        # length of window of event
                    loc = int(ceil((pk_rel/window)*nbins))*20   # phase bin location of amplitude pk
                    
                    # Remove peaks before phase event start / after phase event end
                    if loc<0:
                        loc = nan
                        evt_mask[e] = False
                    elif loc>360:
                        loc = nan
                        evt_mask[e] = False
    
                    # Save phase bin location
                    pkbin[e] = loc
                    
                # Take circular mean of phase locations
                pkbin_rad = asarray([radians(x) for x in pkbin])
                theta = circ_mean(pkbin_rad) # mean circular direction (radians)
                theta_deg = degrees(theta) # mean circular direction (degrees)
                
                # Fix negative angles
                if theta_deg < 0:
                    theta_deg += 360
                
                # Calculate Mean Vector Length (MVL)    
                rad = circ_r(pkbin_rad)
                
                # Rayleigh test for non-uniformity of circular data
                Ray_z, Ray_pv = circ_rayleigh(pkbin_rad)
                Ray_pv = round(Ray_pv,5)
                
    
                ## CALCULATE COUPLING STRENGTH 
                print('')
                print('')
                print('Extracting coupling strength.')
                segs = fetch(dataset, annot, cat=(0,0,0,0), evt_type=[target], cycle=None, 
                                 chan_full=[chan_full], reject_epoch=True, buffer=buffer,
                                 reject_artf = ['Artefact', 'Arou', 'Arousal'])
                segs.read_data([ch], ref_chan, grp_name=grp_name)
                segs_cp = [b for a, b in zip(evt_mask, segs) if a] # mask events (amp peak within phase start&end)
                
                # Initialise variables
                ms = int(ceil(len(segs_cp)/50))
                miraw = zeros((len(events),1))
                longamp = zeros((ms,50),dtype=object) # initialise (blocked) ampltidue series
                longpha = zeros((ms,50),dtype=object) # initialise (blocked) phase series 
                ampbin = zeros((len(segs_cp), nbins)) # initialise mean amplitudes
                z=0 #(progress param)
                for s, seg in enumerate(segs_cp):
                    
                    # Print out progress
                    if progress:
                        z +=1
                        j = z/len(segs_cp)
                        sys.stdout.write('\r')
                        sys.stdout.write(f"Progress: [{'=' * int(50 * j):{50}s}] {int(100 * j)}%")
                        sys.stdout.flush() 
                    
                    # Extract data
                    data = seg['data']
                    
                    # Check polarity of recording
                    if isinstance(polar, list):
                        polarity = polar[i]
                    else:
                        polarity = polar
                    if polarity == 'opposite':
                        data()[0][0] = data()[0][0]*-1   
 
                    
                    # Apply filtering (if necessary)
                    if filter_opts['notch']:
                        selectchans = list(data.chan[0])
                        data.data[0] = notch_mne(data, oREF=filter_opts['oREF'], 
                                                    channel=selectchans, 
                                                    freq=filter_opts['notch_freq'],
                                                    rename=filter_opts['chan_rename'],
                                                    renames=filter_opts['renames'])
                    
                    if filter_opts['notch_harmonics']: 
                        selectchans = list(data.chan[0])
                        data.data[0] = notch_mne2(data, oREF=filter_opts['oREF'], 
                                                  channel=selectchans,
                                                  rename=filter_opts['chan_rename'],
                                                  renames=filter_opts['renames'])
                    
                    if filter_opts['laplacian']:
                        data = laplacian_mne(data, oREF=filter_opts['oREF'], channel=ch, 
                                             ref_chan=ref_chan, 
                                             laplacian_rename=filter_opts['laplacian_rename'], 
                                             renames=filter_opts['renames'])
                        dat = data[0]
                    else:
                        dat = data()[0][0]
                    
                    
                    # Filter data in phase and frequency
                    pha = squeeze(pac.filter(s_freq, dat, ftype='phase'))
                    amp = squeeze(pac.filter(s_freq, dat, ftype='amplitude'))
                    
                    # Put data in blocks (for surrogate testing)
                    longpha[s//50, s%50] = pha
                    longamp[s//50, s%50] = amp
                    
                    # Put data in long format (for preferred phase)
                    ampbin[s, :] = mean_amp(pha, amp, nbins=nbins)
                    
                # If number of events not divisible by block length
                # pad incomplete final block with randomly resampled events
                rem = len(segs_cp) % 50
                if rem > 0:
                    pads = 50 - rem
                    for pad in range(pads):
                        ran = random.randint(0,rem)
                        longpha[-1,rem+pad] = longpha[-1,ran]
                        longamp[-1,rem+pad] = longamp[-1,ran]
                
                # Calculate coupling strength
                mi_r = zeros((longamp.shape[0],1)) # initialise array for raw mi
                mi = zeros((longamp.shape[0],1)) # initialise array for norm mi
                mi_pv = zeros((longamp.shape[0],1)) # initialise array for mi signif.
                for row in range(longamp.shape[0]): 
                    amp = zeros((1))   
                    pha = zeros((1)) 
                    for col in range(longamp.shape[1]):
                        pha = concatenate((pha,longpha[row,col]))
                        amp = concatenate((amp,longamp[row,col]))
                    pha = reshape(pha,(1,1,len(pha)))
                    amp = reshape(amp,(1,1,len(amp)))
                    mi[row] = pac.fit(pha, amp, n_perm=400,random_state=5,
                                 verbose=False)[0][0]
                    mi_r[row] = pac.pac[0][0][0]
                    mi_pv[row] = pac.infer_pvalues(p=0.95, mcp='fdr')[0][0]
                mi_raw = nanmean(mi_r)
                mi_norm = nanmean(mi)
                mi_sig = median(mi_pv)
                
                # Correlation between mean amplitudes and phase-giving sine wave
                # Take mean across all segments/events
                ampbin = ampbin / ampbin.sum(-1, keepdims=True) # normalise amplitude
                ampbin = ampbin.squeeze()
                ma = nanmean(ampbin, axis=0)
                sine = sin(linspace(-pi, pi, nbins))
                sine = interp(sine, (sine.min(), sine.max()), (ma.min(), ma.max()))
                rho, rho_pv = circ_corrcc(ma, sine)
                
                # Add to output dataframe
                out[0,0] = mi_raw
                out[0,1] = mi_norm
                out[0,2] = mi_sig
                out[0,3] = theta
                out[0,4] = theta_deg 
                out[0,5] = rad
                out[0,6] = rho
                out[0,7] = rho_pv
                out[0,8] = Ray_z
                out[0,9] = Ray_pv
    
                # Save cfc metrics to group dataframe
                d = DataFrame(out)
                d.columns = ['mi_raw','mi_norm','mi_sig','pp_rad','ppdegrees','mvl',
                             'rho','pval','ray_z','ray_pv']
                stagename = ''.join(stage) 
                d.to_csv(path_or_buf=out_dir + '/' + p + '/' + vis + '/' + 
                         p + '_' + vis + '_' + ch + '_' + stagename + '_' + band_pairs + 
                         '_cfc_params.csv', sep=',')
                
                # Save binned amplitudes to pickle file
                with open(out_dir + '/' + p + '/' + vis + '/' + 
                         p + '_' + vis + '_' + ch + '_' + stagename + '_' + band_pairs + 
                          '_mean_amps.p', 'wb') as f:
                     dump(ampbin, f)
    
    print('The function pac_it_joint completed without error.')
    return
                            
def cfc_grouplevel(in_dir, out_dir, band_pairs, part, visit, chan, stage, cat,
                   cycle_idx):                        
             
    '''
    This script combines the output from the function pac_it, and formats it
    in a group-level dataframe for statistical analyses.
    The outputs provided by this script will be, for each visit and EEG channel:
        1. A csv array of size:
            i. [#subjects x #phasebins] (if cycle_idx is None), or;
            ii.[#subjects x #sleep cycles x #phasebins] (if cycle_idx is a list)
        2. A csv dataframe with the PAC metrics selected in the analyses from 
            pac_it.  
        
    '''   
    
    # Make output directory
    if not path.exists(out_dir):
        mkdir(out_dir)
    
    ## BIDS CHECKING
    # Check input participants
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(in_dir)
            part = [ p for p in part if not '.' in p]
    else:
        print('')
        print(colored('ERROR |', 'red', attrs=['bold']),
              colored("'part' must either be an array of subject ids or = 'all'" ,
                      'cyan', attrs=['bold']))
        print('')
        
    # Check input visits
    if isinstance(visit, list):
        None
    elif visit == 'all':
        lenvis = set([len(next(walk(in_dir + x))[1]) for x in part])
        if len(lenvis) > 1:
            print(colored('WARNING |', 'yellow', attrs=['bold']),
                  colored('number of visits are not the same for all subjects.',
                          'white', attrs=['bold']))
            print('')
        visit = list(set([y for x in part for y in listdir(in_dir + x)  if '.' not in y]))
    else:
        print('')
        print(colored('ERROR |', 'red', attrs=['bold']),
              colored("'visit' must either be an array of subject ids or = 'visit'" ,
                      'cyan', attrs=['bold']))
        print('')
            
    # Create output dataframe
    if cycle_idx is not None:
        all_ampbin = zeros((7, len(part), 6), dtype='object')
    else:
        all_ampbin = zeros((7, len(part)), dtype='object')
    
    # Check for stage setup
    if cat[1] == 1:
        stage = [''.join(stage) ]

    for st, stagename in enumerate(stage): # Loop through stages
        for k, ch in enumerate(chan):      # Loop through channels
            print('')
            print(f'CHANNEL {ch}')
            for j, vis in enumerate(visit): 
                z=0
                index=[]
                part.sort()
                for i, p in enumerate(part):    # Loop through participants
                    index.append(p)
                    if not path.exists(in_dir + '/' + p + '/' + vis + '/'):
                        print(colored('WARNING |', 'yellow', attrs=['bold']),
                              colored(f'input folder missing for Subject {p}, visit {vis}, skipping..',
                                      'white', attrs=['bold']))
                        continue
                    else:
                        
                        # MEAN AMPLITUDES
                        # Define pickle files for mean amplitudes
                        p_files = [s for s in listdir(in_dir + '/' + p + '/' + vis) 
                                   if band_pairs in s if '.p' in s] 
                        p_files = [s for s in p_files if ch in s]
                        p_files = [s for s in p_files if '_'+stagename+'_' in s]
                        
                        # Open files containing mean amplitudes (if existing)
                        if len(p_files) == 0:
                            print(colored('WARNING |', 'yellow', attrs=['bold']),
                                  colored(f'mean amplitudes file does not exist for Subject {p}, visit {vis}, stage {stagename}, channel {ch}. Skipping..',
                                          'white', attrs=['bold']))
                        elif len(p_files)>1:
                            print(colored('WARNING |', 'yellow', attrs=['bold']),
                                  colored(f'multiple mean amplitudes files exist for Subject {p}, visit {vis}, stage {stagename}, channel {ch}. Skipping..',
                                          'white', attrs=['bold']))
                            print('')
                        else:
                            print(f'Extracting mean amps for ... Subject {p}, visit {vis}')
                            ab_file = in_dir + '/' + p + '/' + vis + '/' + p_files[0]
                            with open(ab_file, 'rb') as f:
                                ampbin = load(f)
                        
                            # Average & normalise mean amplitudes for all segments across the night
                            if cycle_idx is not None:
                                for l in range(0, size(ampbin,0)):
                                    if size(ampbin,0) > 1:
                                        all_ampbin[j, i, l] = nanmean(ampbin[l] / ampbin[l].sum(-1, keepdims=True),
                                                                   axis=0)
                                    else:
                                        all_ampbin[j, i, l] = ampbin[l][0]
                            else:
                                if size(ampbin,0) > 1:
                                        all_ampbin[j, i] = nanmean(ampbin / ampbin.sum(-1, keepdims=True), axis=0)
                                else:
                                    all_ampbin[j, i] = mean(ampbin[0] / ampbin[0].sum(-1, keepdims=True), axis=0)
                        
                        # Define csv files for PAC parameters
                        c_files = [s for s in listdir(in_dir + '/' + p + '/' + vis) 
                                   if band_pairs in s if '.csv' in s]
                        c_files = [s for s in c_files if '_'+stagename+'_' in s] 
                        c_files = [s for s in c_files if ch in s]
                        
                        # Open files containing cfc parameters (if existing)
                        if len(c_files) == 0:
                            print(colored('WARNING |', 'yellow', attrs=['bold']),
                                  colored(f'PAC csv file does not exist for Subject {p}, visit {vis}, stage {stagename}, channel {ch}. Skipping..',
                                          'white', attrs=['bold']))
                            print('')
                        elif len(c_files)>1:
                            print(colored('WARNING |', 'yellow', attrs=['bold']),
                                  colored(f'multiple PAC csv files exist for Subject {p}, visit {vis}, stage {stagename}, channel {ch}. Skipping..',
                                          'white', attrs=['bold']))
                            print('')
                        else:
                            print(f'Extracting PAC params for... Subject {p}, visit {vis}')
                            
                            if z == 0:
                                df = read_csv(in_dir + '/' + p + '/' + vis + '/' + c_files[0],
                                              index_col = 0)
                                df.index = [p]
                                df['chan'] = ch
                                z+=1
                            else:
                                dfload = read_csv(in_dir + '/' + p + '/' + vis + '/' + c_files[0],
                                              index_col = 0)
                                dfload.index = [p]
                                dfload['chan'] = ch
                                df = concat([df, dfload])
                                z+=1
                    
                # Rearrange columns of dataframe   
                try:
                    df = df[[df.columns[-1]] + df.columns[0:-1].tolist()]
    
                    # Save PAC parameters for all participants & visits to file
                    df.to_csv(path_or_buf=out_dir + '/' + ch + '_' + stagename 
                                      + '_visit_' + vis +  '_' + band_pairs + '_cfc_params.csv', 
                                      sep=',')
                
                    
                    # Save mean amplitudes for all participants & visits to file  
                    vis_ampbin = [all_ampbin[j,x] for x in range(0,size(all_ampbin,1))] 
                    vis_ampbin = DataFrame(vis_ampbin, index=index)
                    vis_ampbin.to_csv(path_or_buf=out_dir + '/' + ch + '_' + stagename 
                                      + '_visit_' + vis +  '_' + band_pairs + '_mean_amps.csv', 
                                      sep=',')
                    print('The function cfc_grouplevel completed without error.')  
                    print('')
                except:
                    print(colored('ERROR |', 'red', attrs=['bold']),
                          colored(f'could not find correct PAC csv files.',
                                  'white', attrs=['bold']),
                          colored("1. Check cat variable. 2. Check names of channels, stages..",
                                  'white', attrs=['bold']))
                    print('')
        
                                            
                        

def generate_adap_bands(peaks,width,chan):
    
    '''
    Generates adapted bands of 2 x width from a file containing spectral peaks,
    for the specified channels
    '''
    peaks1 = read_csv(peaks, delimiter=',',index_col=0)
    peaks2 = DataFrame(nan, index=peaks1.index, columns=peaks1.columns)

    
    for c,ch in enumerate(chan):
        peaks2[ch] =  [(x - 2.0, x + 2.0) for x in peaks1[ch]] 
         

    return peaks2    



def watson_williams(in_dir, out_dir, band_pairs, chan, cycle_idx, stage, nbins,
                    test = 'within', comps = [('all','V1'), ('all','V2')]):
    
    '''
    This script conducts a Watson-Williams test between two time-points (within)
    or between 2 groups (between)
            
    '''  
    
    if len(comps)>2:
        print('')
        print('Please only provide 2 comparisons at a time in comps.')
        
    else:
        
        # Setup output directory
        if path.exists(out_dir):
                print(out_dir + " already exists")
        else:
            mkdir(out_dir)
            
        
        # Check if band_pairs is a list
        if isinstance(band_pairs,str):
            band_pairs = [band_pairs]
            
        # Create output file
        dset = zeros((len(chan),len(band_pairs)*2))

        # Set vecbin
        width = 2 * pi / nbins
        vecbin = zeros(nbins)
        for i in range(nbins):
            vecbin[i] = i * width + width / 2
        
        # Loop through channels
        for k, ch in enumerate(chan):
            for b,bp in enumerate(band_pairs):
                
                print('')
                print(f'CHANNEL: {ch}')
                print(f'BAND PAIR: {bp}')
                print('')
                
                # Create output filename    
                stagename = ''.join(stage)
                partstr = ['_'.join(x) for x in comps]
                comparisons = [partstr[0], partstr[1]]
                bands = '_'.join(band_pairs)
                outname = '_vs_'.join([x for x in comparisons])
                filename = f'{stagename}_{bands}_{outname}'
                
                data_m = []
                # Loop through comparisons
                for c,(part,visit) in enumerate(comps):
                
                    # Set list of input participants & visits
                    if isinstance(part, list):
                        None
                    elif part == 'all':
                            part = listdir(in_dir)
                            part = [ p for p in part if not '.' in p]
                    else:
                        print('')
                        print("ERROR: comps must either contain a list of subject ids or = 'all' ")
                        print('')
                    part.sort()
                    for i, p in enumerate(part):
                        if visit == 'all':
                            visit = listdir(in_dir + '/' + p)
                            visit = [x for x in visit if not'.' in x]  
                    if isinstance(visit,str):
                        visit = [visit]
                    visit.sort()    
                    # Define output object   
                    datab = zeros((len(part),len(visit),nbins))
                    
                    # Loop through participants & visits
                    for i, p in enumerate(part):
                        for j, vis in enumerate(visit): 
                            if not path.exists(in_dir + '/' + p + '/' + vis + '/'):
                                print(f'WARNING: input folder missing for Subject {p}, visit {vis}, skipping..')
                                continue
                            else:
                                p_files = [s for s in listdir(in_dir + '/' + p + '/' + vis) if 
                                           bp in s if '.p' in s]
                                p_files = [s for s in p_files if stagename in s] 
                                p_files = [s for s in p_files if ch in s]
                                
                                # Open files containing mean amplitudes (if existing)
                                if len(p_files) == 0:
                                    print(f'WARNING: mean amplitudes file does not exist for Subject {p}, visit {vis}, stage {stagename}, channel {ch} - check this. Skipping..')
                                elif len(p_files) >1:
                                    print(f'WARNING: multiple mean amplitudes files exist for Subject {p}, visit {vis}, stage {stagename}, channel {ch} - check this. Skipping..')
                                else:
                                    print(f'Extracting... Subject {p}, visit {vis}')
                                    ab_file = in_dir + '/' + p + '/' + vis + '/' + p_files[0]
                                    with open(ab_file, 'rb') as f:
                                        ab = load(f)
                                    
                                    # Calculate mean amplitudes across night per subject
                                    ab = nanmean(ab, axis=0)
    
                                    # Caculate z-score for binned data
                                    datab[i,j,:] = zscore(ab)
                                 
                                    
                    # Remove nans from output and take average           
                    databz = array([[[x if not isnan(x) else 0 for x in dim1] 
                                     for dim1 in dim2] for dim2 in datab])
                    data = mean(databz,axis=1)
                    data_m.append(array([circ_mean(vecbin, data[x, :] * 1000) for x in 
                                  range(data.shape[0])]))
                    
                # Create array of data
    
                if test == 'within':
                    if len(data_m[0]) == len(data_m[1]):
                        # Run permutation testing
                        print('')
                        print("Running 10,000 permutations... ")
                        F = zeros((10000))
                        P = zeros((10000))
                        warnings = True
                        for pm in range(0,10000):
                            perm = random.choice(a=[False, True], size=(len(data_m[0])))
                            da = copy.deepcopy(data_m[0])
                            db = copy.deepcopy(data_m[1])
                            if pm>0:
                                da[perm] = data_m[1][perm]
                                db[perm] = data_m[0][perm]
                                warnings = False
                                
                            F[pm], P[pm] = circ_wwtest(da, db, ones(da.shape), 
                                                       ones(db.shape), warnings)
                        dset[k,b*2] = F[0]
                        dset[k,(b*2)+1] = sum(F>F[0])/len(F)
                    else:
                        print("For within-subjects comparisons, the number of subjects in each condition need to match... ")
                elif test == 'between':
                    da = copy.deepcopy(data_m[0])
                    db = copy.deepcopy(data_m[1])
                    F, P = circ_wwtest(da, db, ones(da.shape), ones(db.shape))
                    dset[k,b*2] = F
                    dset[k,(b*2)+1] = P
                else:
                    print("WARNING: test must either be 'between' or 'within' ... ")
            
            
        # Save output to file
        columns = [x+'_'+y for x in band_pairs for y in ['F','p']]
        df = DataFrame(dset, index=chan, columns=columns)
        df.to_csv(r"{out_dir}/watson_williams_{filename}.csv".format(out_dir=out_dir,
                                                           filename=filename))
            
            
            
        print('')
        print("Completed... ")

            

            

