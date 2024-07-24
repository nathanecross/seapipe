# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 15:51:24 2021

@author: Nathan Cross
"""
from datetime import datetime, date

from os import listdir, mkdir, path, walk
from . cfc_func import _allnight_ampbin, circ_wwtest, mean_amp, klentropy
from seapipe.utils.misc import bandpass_mne, laplacian_mne, notch_mne, notch_mne2
from copy import deepcopy
import shutil
from math import degrees, radians
import mne
import matplotlib.pyplot as plt
from numpy import (angle, append, argmax, array, arange, asarray, ceil, concatenate, 
                   empty, histogram, interp, isnan, linspace, log, logical_and, mean, 
                   median, nan, nanmean, ndarray, newaxis, ones, pi, random, repeat, 
                   reshape, roll, save, sin, size, squeeze, sqrt, std, sum, tile, where, zeros) 
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
from seapipe.utils.logs import create_logger, create_logger_outfile
from ..utils.load import (load_channels, load_adap_bands, rename_channels, read_inversion, read_manual_peaks)
from ..utils.misc import remove_duplicate_evts


def pac_method(method, surrogate, correction, list_methods=False):
    
    ''' Formats the method and corrections to be applied through Tensorpac:
        https://etiennecmb.github.io/tensorpac/auto_examples/index.html#tutorials
    '''
    # Calculate Coupling Strength (idpac)
    methods = {1: 'Mean Vector Length (MVL) [Canolty et al. 2006 (Science)]',
               2 : 'Modulation Index (MI) [Tort 2010 (J Neurophys.)]',
               3 : 'Heights Ratio (HR) [Lakatos 2005 (J Neurophys.)]',
               4 : 'ndPAC [Ozkurt 2012 (IEEE)]',
               5 : 'Phase-Locking Value (PLV) [Penny 2008 (J. Neuro. Meth.), Lachaux 1999 (HBM)]',
               6 : 'Gaussian Copula PAC (GCPAC) `Ince 2017 (HBM)`'}
    surrogates = {0 :' No surrogates', 
                  1 : 'Swap phase / amplitude across trials [Tort 2010 (J Neurophys.)]',
                  2 : 'Swap amplitude time blocks [Bahramisharif 2013 (J. Neurosci.) ]',
                  3 : 'Time lag [Canolty et al. 2006 (Science)]'}
    corrections = {0 : 'No normalization',
                   1 : 'Substract the mean of surrogates',
                   2 : 'Divide by the mean of surrogates',
                   3 : 'Substract then divide by the mean of surrogates',
                   4 : 'Z-score'}
    
    if list_methods:
        idpac = [methods] + [surrogates] + [corrections]
    else:
        meth = [x for x in methods if method in methods[x]][0]
        surr = [x for x in surrogates if surrogate in surrogates[x]][0]
        corr = [x for x in corrections if correction in corrections[x]][0]
        idpac = tuple((meth,surr,corr))
    
    return idpac


class pacats:

    def __init__(self, rec_dir, xml_dir, out_dir, log_dir, chan, ref_chan, 
                 grp_name, stage, rater = None, subs = 'all', 
                 sessions = 'all', reject_artf = ['Artefact', 'Arou', 'Arousal'], 
                 tracking = None):
        
        self.rec_dir = rec_dir
        self.xml_dir = xml_dir
        self.out_dir = out_dir
        self.log_dir = log_dir
        
        self.chan = chan
        self.ref_chan = ref_chan
        self.grp_name = grp_name
        self.stage = stage
        self.rater = rater
        self.reject = reject_artf
        
        self.subs = subs
        self.sessions = sessions
        
        if tracking == None:
            tracking = {}
        self.tracking = tracking


    def pac_it(self, cycle_idx, cat, nbins, filter_opts, 
                 epoch_opts, frequency_opts, event_opts, filetype = '.edf',
                 idpac = (2,3,4), evt_type = None, min_dur = 1, 
                 adap_bands_phase = 'Fixed', frequency_phase = (0.5,1.25), 
                 adap_bands_amplitude = 'Fixed', frequency_amplitude = (11,16), 
                 peaks = None, adap_bw = 4, invert = False, progress=True,
                 outfile = 'event_coupling_log.txt'):

        '''
        P.A.C.A.T.S
        
        Phase Amplitude Coupling Across Time Series
        
        This script runs Phase Amplitude Coupling analyses on continuous sleep 
        EEG data (unlike O.C.T.O.P.U.S. which is for events). 
        The method for calculating PAC is set by the parameter <idpac>. 
        For more information on the available methods, refer to the documentation of 
        tensorpac (https://etiennecmb.github.io/tensorpac/index.html) or the article
        (Combrisson et al. 2020, PLoS Comp Bio: https://doi.org/10.1371/journal.pcbi.1008302)
        
        The script does the following:
            1. Extracts the continuous EEG signal segments defined by <cat> and 
                if the segment is too short ± a buffer on either side of length 
                (in sec) specified by <buffer>.
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
        # Get method descriptions
        pac_list = pac_method(0, 0, 0, list_methods=True)
        methods = pac_list[0]
        surrogates = pac_list[1]
        corrections = pac_list[2]
    
        ### 0.a Set up logging
        tracking = self.tracking
        flag = 0
        if outfile == True:
            today = date.today().strftime("%Y%m%d")
            now = datetime.now().strftime("%H:%M:%S")
            logfile = f'{self.log_dir}/event_coupling_{today}_log.txt'
            logger = create_logger_outfile(logfile=logfile, name='Phase-amplitude coupling')
            logger.info('')
            logger.info(f"-------------- New call of 'Phase Amplitude Coupling' evoked at {now} --------------")
        elif outfile:
            logfile = f'{self.log_dir}/{outfile}'
            logger = create_logger_outfile(logfile=logfile, name='Phase-amplitude coupling')
        else:
            logger = create_logger('Phase-amplitude coupling')
        
        logger.info('')
        logger.debug(rf"""Commencing phase-amplitude coupling pipeline... 
                                
                     
                                                                       /~~\
                        ____                                         /'o  |
                     .';;|;;\            _,-;;;\;-_               ,'  _/'|
                    `\_/;;;/;\         /;;\;;;;\;;;,             |     .'
                       `;/;;;|      ,;\;;;|;;;|;;;|;\          ,';;\  |
                       |;;;/;:     |;;;\;/~~~~\;/;;;|        ,;;;;;;.'
                      |;/;;;|     |;;;,'      `\;;/;|      /;\;;;;/
                      `|;;;/;\___/;~\;|         |;;;;;----\;;;|;;/'
                       `;/;;;|;;;|;;;,'         |;;;;|;;;;;|;;|/'
                        `\;;;|;;;/;;,'           `\;/;;;;;;|/~'
                         `\/;;/;;;/               `~------'
                           `~~~~~  

                
                Phase Amplitude Coupling Across Time Series
                (P.A.C.A.T.S)
                
                Method: {methods[idpac[0]]}
                Correction: {surrogates[idpac[1]]}
                Normalisation: {corrections[idpac[2]]}
                                  
                                                    """,)
        ### 0.b. Set up organisation of export
        if cat[0] + cat[1] == 2:
            model = 'whole_night'
            logger.debug('Analysing PAC for the whole night.')
        elif cat[0] + cat[1] == 0:
            model = 'stage*cycle'
            logger.debug('Analysing PAC per stage and cycle separately.')
        elif cat[0] == 0:
            model = 'per_cycle'
            logger.debug('Analysing PAC per cycle separately.')
        elif cat[1] == 0:
            model = 'per_stage'  
            logger.debug('Analysing PAC per stage separately.')
        if 'cycle' in model and cycle_idx == None:
            logger.info('')
            logger.critical("To run cycles separately (i.e. cat[0] = 0), cycle_idx cannot be 'None'")
            return
        
        # Log filtering options
        if filter_opts['notch']:
            logger.debug(f"Applying notch filtering: {filter_opts['notch_freq']} Hz")
        if filter_opts['notch_harmonics']: 
            logger.debug('Applying notch harmonics filtering.')
        if filter_opts['bandpass']:
            logger.debug(f"Applying bandpass filtering: {filter_opts['highpass']} - {filter_opts['lowpass']} Hz")
        if filter_opts['laplacian']:
            logger.debug('Applying Laplacian filtering.')
        
        ### 1. First we check the directories
        # a. Check for output folder, if doesn't exist, create
        if path.exists(self.out_dir):
                logger.debug("Output directory: " + self.out_dir + " exists")
        else:
            mkdir(self.out_dir)
        
        # b. Check input list
        subs = self.subs
        if isinstance(subs, list):
            None
        elif subs == 'all':
                subs = listdir(self.rec_dir)
                subs = [p for p in subs if not '.' in p]
        else:
            logger.error("'subs' must either be an array of subject ids or = 'all' ")   
        
        ### 2. Begin loop through dataset
        # a. Begin loop through participants
        subs.sort()
        for i, sub in enumerate(subs):
            tracking[f'{sub}'] = {}
            # b. Begin loop through sessions
            sessions = self.sessions
            if sessions == 'all':
                sessions = listdir(self.rec_dir + '/' + sub)
                sessions = [x for x in sessions if not '.' in x]   
            
            for v, ses in enumerate(sessions):
                logger.info('')
                logger.debug(f'Commencing {sub}, {ses}')
                tracking[f'{sub}'][f'{ses}'] = {'pac':{}} 
                
                
                ## c. Load recording
                rdir = self.rec_dir + '/' + sub + '/' + ses + '/eeg/'
                try:
                    edf_file = [x for x in listdir(rdir) if x.endswith(filetype)]
                    dset = Dataset(rdir + edf_file[0])
                except:
                    logger.warning(f' No input {filetype} file in {rdir}')
                    break
                
                ## d. Load annotations
                xdir = self.xml_dir + '/' + sub + '/' + ses + '/'
                try:
                    xml_file = [x for x in listdir(xdir) if x.endswith('.xml')]
                    # Copy annotations file before beginning
                    if not path.exists(self.out_dir):
                        mkdir(self.out_dir)
                    if not path.exists(self.out_dir + '/' + sub):
                        mkdir(self.out_dir + '/' + sub)
                    if not path.exists(self.out_dir + '/' + sub + '/' + ses):
                        mkdir(self.out_dir + '/' + sub + '/' + ses)
                    outpath = self.out_dir + '/' + sub + '/' + ses
                    backup_file = (f'{outpath}/{sub}_{ses}_spindle.xml')
                    if not path.exists(backup_file):
                        shutil.copy(xdir + xml_file[0], backup_file)
                    else:
                        logger.debug(f'Using annotations file from: {xdir}')
                except:
                    logger.warning(f' No input annotations file in {xdir}')
                    flag +=1
                    break
                
                # 3.a. Read annotations file
                annot = Annotations(backup_file, rater_name=self.rater)
                
                ## b. Get sleep cycles (if any)
                if cycle_idx is not None:
                    all_cycles = annot.get_cycles()
                    cycle = [all_cycles[y - 1] for y in cycle_idx if y <= len(all_cycles)]
                else:
                    cycle = None
                
                ## c. Channel setup 
                pflag = deepcopy(flag)
                flag, chanset = load_channels(sub, ses, self.chan, self.ref_chan,
                                              flag, logger)
                if flag - pflag > 0:
                    logger.warning(f'Skipping {sub}, {ses}...')
                    flag +=1
                    break
                newchans = rename_channels(sub, ses, self.chan, logger)
                
                # 4.a. Loop through channels
                for c, ch in enumerate(chanset):
                    
                    # b. Rename channel for output file (if required)
                    if newchans:
                        fnamechan = newchans[ch]
                    else:
                        fnamechan = ch

                    if ch == '_REF':
                        filter_opts['oREF'] = newchans[ch]
                    else:
                        filter_opts['oREF'] = None
                        
                    # c. Set frequency bands for:
                    # Phase
                    if adap_bands_phase == 'Fixed':
                        f_pha = frequency_phase   
                    elif adap_bands_phase == 'Manual':
                        f_pha = read_manual_peaks(sub, ses, peaks, ch, 
                                                 adap_bw, logger)
                    elif adap_bands_phase == 'Auto':
                        stagename = '-'.join(self.stage)
                        band_limits = f'{self.frequency[0]}-{self.frequency[1]}Hz'
                        f_pha = load_adap_bands(self.tracking['fooof'], sub, ses,
                                               fnamechan, stagename, band_limits, 
                                               adap_bw, logger)
                    if not f_pha:
                        logger.warning('Will use fixed frequency bands for PHASE instead.')
                        f_pha = frequency_phase
                    if not chanset[ch]:
                        logchan = ['(no re-refrencing)']
                    else:
                        logchan = chanset[ch]
                    logger.debug(f"Using PHASE frequency band: {round(f_pha[0],2)}-{round(f_pha[1],2)} Hz for {sub}, {ses}, {str(ch)}:{'-'.join(logchan)}")    
                    
                    # Amplitude
                    if adap_bands_amplitude == 'Fixed':
                        f_amp = frequency_amplitude   
                    elif adap_bands_amplitude == 'Manual':
                        f_amp = read_manual_peaks(sub, ses, peaks, ch, 
                                                 adap_bw, logger)
                    elif adap_bands_amplitude == 'Auto':
                        stagename = '-'.join(self.stage)
                        band_limits = f'{self.frequency[0]}-{self.frequency[1]}Hz'
                        f_amp = load_adap_bands(self.tracking['fooof'], sub, ses,
                                               fnamechan, stagename, band_limits, 
                                               adap_bw, logger)
                    if not f_amp:
                        logger.warning('Will use fixed frequency bands for PHASE instead.')
                        f_amp = frequency_amplitude
                    if not chanset[ch]:
                        logchan = ['(no re-refrencing)']
                    else:
                        logchan = chanset[ch]
                    logger.debug(f"Using AMPLITUDE frequency band: {round(f_amp[0],2)}-{round(f_amp[1],2)} Hz for {sub}, {ses}, {str(ch)}:{'-'.join(logchan)}")    
                    
                    # d. Check if channel needs to be inverted for detection
                    if type(invert) == type(DataFrame()):
                        inversion = read_inversion(sub, ses, invert, ch, logger)
                        if not inversion:
                            logger.warning(f"NO inversion will be applied to channel {ch} prior to detection for {sub}, {ses}. To turn off this warning, select `invert = 'False'`")
                        else: 
                            logger.debug(f'Inverting channel {ch} prior to detection for {sub}, {ses}')
                    elif type(invert) == bool:
                        inversion = invert
                    else:
                        logger.critical(f"The argument 'invert' must be set to either: 'True', 'False' or 'None'; but it was set as {invert}.")
                        logger.info('Check documentation for how to set up staging data:')
                        logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
                        logger.info('-' * 10)
                        return
    
                    # 5.a. Fetch data
                    logger.debug(f"Reading EEG data for {sub}, {ses}, {str(ch)}:{'-'.join(logchan)}")
                    if not evt_type == None and not isinstance(evt_type, list):
                        evt_type = [evt_type]
                    segments = fetch(dset, annot, cat = cat, evt_type = evt_type, 
                                     stage = self.stage,  cycle = cycle, 
                                     buffer = event_opts['buffer'])
                    # segments = fetch(dset, annot, cat = cat, evt_type = evt_type, 
                    #                  stage = self.stage, cycle=cycle,   
                    #                  epoch = epoch_opts['epoch'], 
                    #                  epoch_dur = epoch_opts['epoch_dur'], 
                    #                  epoch_overlap = epoch_opts['epoch_overlap'], 
                    #                  epoch_step = epoch_opts['epoch_step'], 
                    #                  reject_epoch = epoch_opts['reject_epoch'], 
                    #                  reject_artf = epoch_opts['reject_artf'],
                    #                  min_dur = epoch_opts['min_dur'])
                    if len(segments)==0:
                        logger.warning(f"No valid data found for {sub}, {ses}, {self.stage}, Cycles:{cycle}.")
                        flag +=1
                        break
                    
                    # 5.b. Read data
                    if filter_opts['laplacian']:
                        try:
                            segments.read_data(filter_opts['lapchan'], chanset[ch]) 
                            laplace_flag = True
                        except:
                            logger.error(f"Channels listed in filter_opts['lapchan']: {filter_opts['lapchan']} are not found in recording for {sub}, {ses}.")
                            logger.warning("Laplacian filtering will NOT be run for {sub}, {ses}, {ch}. Check parameters under: filter_opts['lapchan']")
                            segments.read_data(ch, chanset[ch])
                            laplace_flag = False
                            flag += 1
                    else:
                        segments.read_data(ch, chanset[ch])
                    
                    
                    # 6.a. Define PAC object
                    pac = Pac(idpac = idpac, f_pha = f_pha, f_amp = f_amp, 
                              dcomplex = filter_opts['dcomplex'], 
                              cycle = filter_opts['filtcycle'], 
                              width = filter_opts['width'], 
                              n_bins = nbins,
                              verbose='ERROR')
                    
                    # b. Create blocks
                    ampbin = zeros((len(segments), nbins))
                    ms = int(ceil(len(segments)/50))
                    longamp = zeros((ms,50),dtype=object) # initialise (blocked) ampltidue series
                    longpha = zeros((ms,50),dtype=object) # initialise (blocked) phase series 
                    
                    # 7.a. Divide segments based on concatenation
                    nsegs=[]
                    if model == 'whole_night':
                        nsegs = [s for s in segments]
                    elif model == 'stage*cycle':
                        for st in self.stage:
                            for cy in cycle_idx:
                                segs = [s for s in segments if st in s['stage'] if cy in s['cycle']]
                                nsegs.append(segs)
                    elif model == 'per_cycle':
                        for cy in cycle_idx:
                            segs = [s for s in segments if cy in s['cycle']]
                            nsegs.append(segs)
                    elif model == 'per_stage':
                        for st in self.stage:
                            segs = [s for s in segments if st in s['stage']]
                            nsegs.append(segs)
                    
                    # b. Loop over segments and apply filtering (if required)
                    z=0
                    logger.info('')
                    for sg in range(len(nsegs)):
                        
                        # Print out progress
                        if progress:
                            z +=1
                            j = z/len(segments)
                            sys.stdout.write('\r')
                            sys.stdout.write(f"                      Progress: [{'»' * int(50 * j):{50}s}] {int(100 * j)}%")
                            sys.stdout.flush()
                        
                        seg = nsegs[sg]
                        out = dict(seg)
                        data = seg['data']
                        timeline = data.axis['time'][0]
                        out['start'] = timeline[0]
                        out['end'] = timeline[-1]
                        out['duration'] = len(timeline) / data.s_freq
                        if filter_opts['laplacian']:
                            selectchans = filter_opts['lapchan']
                        else:
                            selectchans = ch
                        
                        # b. Notch filters
                        if filter_opts['notch']:
                            data.data[0] = notch_mne(data, oREF=filter_opts['oREF'], 
                                                        channel=selectchans, 
                                                        freq=filter_opts['notch_freq'],
                                                        rename=filter_opts['laplacian_rename'],
                                                        renames=filter_opts['renames'],
                                                        montage=filter_opts['montage'])
                            
                        if filter_opts['notch_harmonics']: 
                            data.data[0] = notch_mne2(data, oREF=filter_opts['oREF'], 
                                                      channel=selectchans, 
                                                      rename=filter_opts['laplacian_rename'],
                                                      renames=filter_opts['renames'],
                                                      montage=filter_opts['montage'])    
                        
                        # c. Bandpass filters
                        if filter_opts['bandpass']:
                            data.data[0] = bandpass_mne(data, oREF=filter_opts['oREF'], 
                                                      channel=selectchans,
                                                      highpass=filter_opts['highpass'], 
                                                      lowpass=filter_opts['lowpass'], 
                                                      rename=filter_opts['laplacian_rename'],
                                                      renames=filter_opts['renames'],
                                                      montage=filter_opts['montage'])
                        
                        # d. Laplacian transform
                        if filter_opts['laplacian'] and laplace_flag:
                            data.data[0] = laplacian_mne(data, 
                                                 filter_opts['oREF'], 
                                                 channel=selectchans, 
                                                 ref_chan=chanset[ch], 
                                                 laplacian_rename=filter_opts['laplacian_rename'], 
                                                 renames=filter_opts['renames'],
                                                 montage=filter_opts['montage'])
                            data.axis['chan'][0] = asarray([x for x in chanset])
                            selectchans = ch
                            dat = data[0]
                        else:
                            dat = data()[0][0]
    
    
                       
                        # e. Fix polarity of recording
                        if inversion:
                            dat = dat*-1 
       
                        # f. Obtain phase signal
                        pha = squeeze(pac.filter(data.s_freq, dat, ftype='phase'))
                        
                        if len(pha.shape)>2:
                            pha = squeeze(pha)
                        
                        # g. obtain amplitude signal
                        amp = squeeze(pac.filter(data.s_freq, dat, ftype='amplitude'))
                        if len(amp.shape)>2:
                            amp = squeeze(amp)
                        
                        # h. extract signal (minus buffer)
                        nbuff = int(event_opts['buffer'] * data.s_freq)
                        minlen = data.s_freq * min_dur
                        if len(pha) >= 2 * nbuff + minlen:
                            pha = pha[nbuff:-nbuff]
                            amp = amp[nbuff:-nbuff]                               
                            
                        # Apply phase correction for hilbert transform
                        #pha = roll(pha, int(pi/2*s_freq), axis=-1)
    
                        # i. put data in blocks (for surrogate testing)
                        longpha[sg//50, sg%50] = pha
                        longamp[sg//50, sg%50] = amp
                        
                        # j. put data in long format (for preferred phase)
                        ampbin[sg, :] = mean_amp(pha, amp, nbins=nbins)
    
                    # 8.a. if number of events not divisible by block length,
                    #    pad incomplete final block with randomly resampled events
                    sys.stdout.write('\r')
                    sys.stdout.flush()
                    
                    rem = len(segments) % 50
                    if rem > 0:
                        pads = 50 - rem
                        for pad in range(pads):
                            ran = random.randint(0,rem)
                            longpha[-1,rem+pad] = longpha[-1,ran]
                            longamp[-1,rem+pad] = longamp[-1,ran]
                    
                    # b. Calculate Coupling Strength
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

                    ## c. Calculate preferred phase
                    ampbin = ampbin / ampbin.sum(-1, keepdims=True) # normalise amplitude
                    ampbin = ampbin.squeeze()
                    ampbin = ampbin[~isnan(ampbin[:,0]),:] # remove nan trials
                    ab = ampbin
                    
                    # d. Create bins for preferred phase
                    vecbin = zeros(nbins)
                    width = 2 * pi / nbins
                    for n in range(nbins):
                        vecbin[n] = n * width + width / 2  
                    
                    # e. Calculate mean direction (theta) & mean vector length (rad)
                    ab_pk = argmax(ab,axis=1)
                    theta = circ_mean(vecbin,histogram(ab_pk,bins=nbins, 
                                                        range=(0,nbins))[0])
                    theta_deg = degrees(theta)
                    if theta_deg < 0:
                        theta_deg += 360
                    rad = circ_r(vecbin, histogram(ab_pk,bins=nbins)[0], d=width)
                    
                    # f. Take mean across all segments/events
                    ma = nanmean(ab, axis=0)
                    
                    # g. Correlation between mean amplitudes and phase-giving sine wave
                    sine = sin(linspace(-pi, pi, nbins))
                    sine = interp(sine, (sine.min(), sine.max()), (ma.min(), ma.max()))
                    rho, pv1 = circ_corrcc(ma, sine)

                    # h. Rayleigh test for non-uniformity of circular data
                    ppha = vecbin[ab.argmax(axis=-1)]
                    z, pv2 = circ_rayleigh(ppha)
                    pv2 = round(pv2,5)
                    
                    # 9.a Export and save data
                    freqs = f'pha-{f_amp[0]}-{f_amp[1]}Hz_amp-{f_amp[0]}-{f_amp[1]}Hz'
                    if model == 'whole_night':
                        stagename = '-'.join(self.stage)
                        outputfile = '{}/{}_{}_{}_{}_{}_pac_parameters.csv'.format(
                                        outpath,sub,ses,fnamechan,stagename,freqs)
                    elif model == 'stage*cycle':    
                        outputfile = '{}/{}_{}_{}_{}_cycle{}_{}_pac_parameters.csv'.format(
                                      outpath,sub,ses,fnamechan,self.stage[sg],cycle_idx[sg],freqs)
                    elif model == 'per_stage':
                        outputfile = '{}/{}_{}_{}_{}_{}_pac_parameters.csv'.format(
                                      outpath,sub,ses,fnamechan,self.stage[sg],freqs)
                    elif model == 'per_cycle':
                        stagename = '-'.join(self.stage)
                        outputfile = '{}/{}_{}_{}_{}_cycle{}_{}_pac_parameters.csv'.format(
                                      outpath,sub,ses,fnamechan,stagename,cycle_idx[sg],freqs)

                    # b. Save cfc metrics to dataframe
                    d = DataFrame([mean(pac.pac), mean(mi), median(mi_pv), theta, 
                                    theta_deg, rad, rho, pv1, z, pv2])
                    d = d.transpose()
                    d.columns = ['mi','mi_norm','sig','pp_rad','ppdegrees','mvl',
                                  'rho','pval','rayl','pval2']
                    d.to_csv(path_or_buf=outputfile, sep=',')
                    
                    # c. Save binned amplitudes to pickle file
                    outputfile = outputfile.split('_pac_parameters.csv')[0] + '_mean_amps'
                    save(outputfile, ab)
       
        ### 10. Check completion status and print
        if flag == 0:
            logger.debug('Phase-amplitude coupling finished without ERROR.')  
        else:
            logger.warning('Phase-amplitude coupling finished with WARNINGS. See log for details.')
        
        #self.tracking = tracking   ## TO UPDATE - FIX TRACKING
        
        return                             


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

            

            

