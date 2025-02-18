# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 15:51:24 2021

@author: Nathan Cross
"""
from datetime import datetime, date

from os import listdir, mkdir, path, walk
from . cfc_func import _allnight_ampbin, circ_wwtest, mean_amp, klentropy
from . octopus import pac_method
from seapipe.utils.misc import bandpass_mne, laplacian_mne, notch_mne, notch_mne2
from copy import deepcopy
import shutil
from math import degrees, radians
import mne
import matplotlib.pyplot as plt
from numpy import (angle, append, argmax, array, arange, asarray, ceil, concatenate, 
                   empty, histogram, interp, isnan, linspace, log, logical_and, mean, 
                   median, nan, nanmean, nanmedian, ndarray, newaxis, ones, pi, random, repeat, 
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
from ..utils.load import (load_channels, load_adap_bands, rename_channels, 
                          load_sessions, read_inversion, read_manual_peaks)
from ..utils.misc import remove_duplicate_evts



class pacats:

    def __init__(self, rootpath, rec_dir, xml_dir, out_dir, chan, ref_chan, 
                 grp_name, stage, rater = None, subs = 'all', 
                 sessions = 'all', reject_artf = ['Artefact', 'Arou', 'Arousal'], 
                 tracking = None):
        
        self.rootpath = rootpath
        self.rec_dir = rec_dir
        self.xml_dir = xml_dir
        self.out_dir = out_dir
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
                 epoch_opts, frequency_opts, filetype = '.edf',
                 idpac = (2,3,4), min_dur = 1, 
                 adap_bands_phase = 'Fixed', frequency_phase = (0.5,1.25), 
                 adap_bands_amplitude = 'Fixed', frequency_amplitude = (11,16), 
                 peaks = None, adap_bw = 4, invert = False, progress=True,
                 logger = create_logger('Phase-amplitude coupling')):

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
    
        tracking = self.tracking
        flag = 0
        
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
            flag, sessions = load_sessions(sub, self.sessions, self.rec_dir, flag, 
                                     logger, verbose=2) 
            
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
                        f_pha = read_manual_peaks(self.roothpath, sub, ses, ch, 
                                                 adap_bw, logger)
                    elif adap_bands_phase == 'Auto':
                        stagename = '-'.join(self.stage)
                        band_limits = f'{frequency_phase[0]}-{frequency_phase[1]}Hz'
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
                        f_amp = read_manual_peaks(self.roothpath, sub, ses, ch, 
                                                 adap_bw, logger)
                    elif adap_bands_amplitude == 'Auto':
                        stagename = '-'.join(self.stage)
                        band_limits = f'{frequency_amplitude[0]}-{frequency_amplitude[1]}Hz'
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
                    segments = fetch(dset, annot, cat = cat, stage = self.stage,  
                                     cycle = cycle)

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
                    
                    
                    # 6. Define PAC object
                    pac = Pac(idpac = idpac, f_pha = f_pha, f_amp = f_amp, 
                              dcomplex = filter_opts['dcomplex'], 
                              cycle = filter_opts['filtcycle'], 
                              width = filter_opts['width'], 
                              n_bins = nbins,
                              verbose='ERROR')
                    
                    # 7.a. Divide segments based on concatenation
                    nsegs=[]
                    if model == 'whole_night':
                        nsegs = [s for s in segments]
                        seg_label = ['whole_night']
                    elif model == 'stage*cycle':
                        seg_label = []
                        for st in self.stage:
                            for cy in cycle_idx:
                                segs = [s for s in segments if st in s['stage'] if cy in s['cycle']][0]
                                nsegs.append(segs)
                                seg_label.append(f'{st}_cycle{cy}')
                    elif model == 'per_cycle':
                        seg_label = []
                        for cy in cycle_idx:
                            segs = [s for s in segments if cy in s['cycle']][0]
                            nsegs.append(segs)
                            seg_label.append(f'cycle{cy}')
                    elif model == 'per_stage':
                        seg_label = []
                        for st in self.stage:
                            segs = [s for s in segments if st in s['stage']][0]
                            nsegs.append(segs)
                            seg_label.append(f'{st}')
                    
                    # 7.b. Loop over segments
                    logger.info('')
                    for nsg in range(len(nsegs)):
                        seg = nsegs[nsg]
                        logger.debug(f'Analysing {seg_label[nsg]}')
                        
                        # 7.c. Create blocks
                        out = dict(seg)
                        data = seg['data']
                        timeline = data.axis['time'][0]
                        out['start'] = timeline[0]
                        out['end'] = timeline[-1]
                        out['duration'] = len(timeline) / data.s_freq

                        epoch = int(ceil((50/mean(frequency_phase))*data.s_freq)) # epoch length of 50 phase cycles
                        ms = int(ceil( len(timeline) / epoch )) # number of epochs in data
                        longamp = zeros((ms, epoch)) # initialise (blocked) amplitude series
                        longpha = zeros((ms, epoch)) # initialise (blocked) phase series 
                        ampbin = zeros((ms, nbins)) # initialise mean
                        
                        
                        
                        # 7.d Apply filtering (if required)
                        if filter_opts['laplacian']:
                            selectchans = filter_opts['lapchan']
                        else:
                            selectchans = ch
                        
                        # 7.d. Notch filters
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
                        
                        # 7.e. Bandpass filters
                        if filter_opts['bandpass']:
                            data.data[0] = bandpass_mne(data, oREF=filter_opts['oREF'], 
                                                      channel=selectchans,
                                                      highpass=filter_opts['highpass'], 
                                                      lowpass=filter_opts['lowpass'], 
                                                      rename=filter_opts['laplacian_rename'],
                                                      renames=filter_opts['renames'],
                                                      montage=filter_opts['montage'])
                        
                        # 7.f. Laplacian transform
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
                            
                            
                        # 7.g. Fix polarity of recording
                        if inversion:
                            dat = dat*-1 
                        
                        for s in range(ms):
                            
                            datseg = dat[s*epoch:(s+1)*epoch]
                            
                            #  pad incomplete final block with nan
                            if len(datseg) < epoch:
                                datseg = concatenate((datseg, empty(epoch-len(datseg))))
                            
                            # Print out progress
                            if progress:
                                j = s/ms
                                sys.stdout.write('\r')
                                sys.stdout.write(f"                      Progress: [{'»' * int(50 * j):{50}s}] {int(100 * j)}%")
                                sys.stdout.flush()
                            
                            # 7.h. Obtain phase signal
                            pha = squeeze(pac.filter(data.s_freq, datseg, ftype='phase'))
                            
                            if len(pha.shape)>2:
                                pha = squeeze(pha)
                            
                            # 7.i. obtain amplitude signal
                            amp = squeeze(pac.filter(data.s_freq, datseg, ftype='amplitude'))
                            if len(amp.shape)>2:
                                amp = squeeze(amp)
                            
                            # 7.j.  put data in blocks (for surrogate testing)
                            longpha[s] = pha
                            longamp[s] = amp
                            
                            # 7.k. put averaged data in list (for preferred phase)
                            ampbin[s, :] = mean_amp(pha, amp, nbins=nbins)
        
                        # 8.a. if number of events not divisible by block length,
                        #    pad incomplete final block with randomly resampled events
                        sys.stdout.write('\r')
                        sys.stdout.flush()
                        
                        # b. Calculate Coupling Strength
                        mi = zeros((longamp.shape[0],1))
                        mi_pv = zeros((longamp.shape[0],1))
                        for r,row in enumerate(longamp): 
                            amp = row
                            pha = longpha[r]
                            amp = amp[~isnan(amp)]
                            if len(amp) < 1:
                                mi[r] = nan
                                mi_pv[r] = nan
                            else:
                                pha = reshape(pha,(1,1,len(pha)))
                                amp = reshape(amp,(1,1,len(amp)))
                                mi[r] = pac.fit(pha, amp, n_perm=400,random_state=5,
                                              verbose=False)[0][0]
                                mi_pv[r] = pac.infer_pvalues(p=0.95, mcp='fdr')[0][0]
    
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
                        if adap_bands_phase == 'Fixed':
                            phadap = '-fixed'
                        else:
                            phadap = '-adap'
                        if adap_bands_amplitude == 'Fixed':
                            ampadap = '-fixed'
                        else:
                            ampadap = '-adap'    
                        freqs = f'pha-{f_pha[0]}-{f_pha[1]}Hz{phadap}_amp-{f_amp[0]}-{f_amp[1]}Hz{ampadap}'
                        if model == 'whole_night':
                            stagename = '-'.join(self.stage)
                            outputfile = '{}/{}_{}_{}_{}_{}_pac_parameters.csv'.format(
                                            outpath,sub,ses,fnamechan,stagename,freqs)
                        elif model == 'stage*cycle':    
                            outputfile = '{}/{}_{}_{}_{}_cycle{}_{}_pac_parameters.csv'.format(
                                          outpath,sub,ses,fnamechan,self.stage[nsg],cycle_idx[nsg],freqs)
                        elif model == 'per_stage':
                            outputfile = '{}/{}_{}_{}_{}_{}_pac_parameters.csv'.format(
                                          outpath,sub,ses,fnamechan,self.stage[nsg],freqs)
                        elif model == 'per_cycle':
                            stagename = '-'.join(self.stage)
                            outputfile = '{}/{}_{}_{}_{}_cycle{}_{}_pac_parameters.csv'.format(
                                          outpath,sub,ses,fnamechan,stagename,cycle_idx[nsg],freqs)
    
                        # b. Save cfc metrics to dataframe
                        d = DataFrame([nanmean(pac.pac), nanmean(mi), nanmedian(mi_pv), theta, 
                                        theta_deg, rad, rho, pv1, z, pv2])
                        d = d.transpose()
                        d.columns = ['mi_raw','mi_norm','pval1','pp_radians','ppdegrees','mvl',
                                      'rho','pval2','rayl','pval3']
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




            

            

