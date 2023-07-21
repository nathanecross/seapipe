# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 10:29:11 2021

@author: Nathan Cross
"""

from csv import reader
from fooof import FOOOF
from fooof.analysis import get_band_peak_fm
from glob import glob
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from numpy import (asarray, ceil, concatenate, empty, floor, mean, nan, 
                   ones, pi, sqrt, stack, sum, zeros)
from openpyxl import Workbook
from os import listdir, mkdir, path
from pandas import DataFrame, ExcelWriter
from scipy.fftpack import next_fast_len
from wonambi import ChanFreq, Dataset
from wonambi.attr import Annotations
from wonambi.trans import (fetch, frequency, get_descriptives, export_freq, 
                           export_freq_band)
from utils.misc import bandpass_mne, laplacian_mne, notch_mne, notch_mne2

def powerspec_it(rec_dir, xml_dir, out_dir, cat, stage, chan, ref_chan, 
                 general_opts, frequency_opts, rater=None, cycle_idx=None, 
                 part='all',visit='all', filter_opts=None, epoch_opts=None, 
                 event_opts=None, norm=None, norm_opts=None, fooof_it = False,
                 fooof_opts=None):
    
    print(r"""    Calculating power spectrum...
          
              |
              | /\ 
              |/  \
          uV2 |    ^-_
              |       ^-___-__
              |________________
                     (Hz)
          
          """)
        
    if path.exists(out_dir):
            print(out_dir + " already exists")
    else:
        mkdir(out_dir)
    
    
    suffix = general_opts['suffix']
        
    if fooof_it:
        if fooof_opts is None:
            print('Error: Options not set for FOOOF')
        else:
            # Define model parameters
            fm = FOOOF(fooof_opts['peak_width_limits'], fooof_opts['max_n_peaks'], 
                       fooof_opts['min_peak_amplitude'], fooof_opts['peak_threshold'])
        
            def gaussian_integral(a, c):
                """ Returns definite integral of a gaussian function with height a and
                standard deviation c."""
                return sqrt(2) * a * abs(c) * sqrt(pi)
        
        # Prepare Bands
        if fooof_opts['bands_fooof'] is not None:
            bands = fooof_opts['bands_fooof']
        else:
            stp = min(fooof_opts['peak_width_limits'])
            low = int(floor(fooof_opts['freq_range'][0]))
            hi = int(ceil(fooof_opts['freq_range'][1]))
            bands = [(x,x + stp) for x in range(low,hi)]
        
    
    
    ### loop through records
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(rec_dir)
            part = [ p for p in part if not '.' in p]
    else:
        print("ERROR: 'part' must either be an array of subject ids or = 'all' ")  

    if norm:
        if norm not in ('integral', 'baseline'):
            exit('Invalid value for norm: ' + str(norm))
            
    for i, p in enumerate(part):
        # loop through visits
        if visit == 'all':
            visit = listdir(rec_dir + '/' + p)
            visit = [x for x in visit if not '.' in x]    
        
        for v, vis in enumerate(visit):
            ## Define files
            rdir = rec_dir + p + '/' + vis + '/'
            xdir = xml_dir + p + '/' + vis + '/'
            edf_file = [x for x in listdir(rdir) if x.endswith('.edf') or x.endswith('.rec') 
                        or x.endswith('.eeg') or x.endswith('.set')]
            xml_file = [x for x in listdir(xdir) if x.endswith('.xml')]            
            
            ### Define output path
            if not path.exists(out_dir):
                mkdir(out_dir)
            if not path.exists(out_dir + p ):
                mkdir(out_dir + p)
            if not path.exists(out_dir + p + '/' + vis):
                mkdir(out_dir + p + '/' + vis)
            outpath = out_dir + p + '/' + vis + '/'
            
            
            ## Now import data
            dset = Dataset(rdir + edf_file[0])
            annot = Annotations(xdir + xml_file[0], rater_name=rater)
             
            ### get cycles
            if cycle_idx is not None:
                all_cycles = annot.get_cycles()
                cycle = [all_cycles[i - 1] for i in cycle_idx if i <= len(all_cycles)]
            else:
                cycle = None
                
            ### if event channel only, specify event channels
            chan_full = None
            if event_opts['event_chan']:
                chan_full = [i + ' (' + general_opts['chan_grp_name'] + ')' for i in event_opts['event_chan']]
            
            ### select and read data
            print('Reading data for ' + p + ', Visit: ' + vis)
            segments = fetch(dset, annot, cat=cat, evt_type=event_opts['evt_type'], 
                             stage=stage, cycle=cycle, chan_full=chan_full, 
                             epoch=epoch_opts['epoch'], 
                             epoch_dur=epoch_opts['epoch_dur'], 
                             epoch_overlap=epoch_opts['epoch_overlap'], 
                             epoch_step=epoch_opts['epoch_step'], 
                             reject_epoch=general_opts['reject_epoch'], 
                             reject_artf=general_opts['reject_artf'],
                             min_dur=general_opts['min_dur'])
            
            if not segments:
                print('No valid data found. Skipping to next record.')
                continue
            
            # BASELINE NORMALISATION
            if norm == 'baseline':
                norm_seg = fetch(dset, annot, cat=norm_opts['norm_cat'], 
                                 evt_type=norm_opts['norm_evt_type'],
                                 stage=norm_opts['norm_stage'], 
                                 epoch=norm_opts['norm_epoch'])
                
                if not norm_seg:
                    print('No valid normalization data found. '
                          'Skipping to next record.')
                    continue
                
                if filter_opts['laplacian']:
                    norm_seg.read_data(filter_opts['lapchan'], ref_chan) 
                else:
                    norm_seg.read_data(chan, ref_chan)            
                all_nSxx = []
                
                for seg in norm_seg:
                    
                    normdata = seg['data']
                    
                    if filter_opts['laplacian']:
                        normdata.data[0] = laplacian_mne(normdata, 
                                                 filter_opts['oREF'], 
                                                 channel=chan, 
                                                 ref_chan=ref_chan, 
                                                 laplacian_rename=filter_opts['laplacian_rename'], 
                                                 renames=filter_opts['renames'])
                    
                    Sxx = frequency(normdata, output=frequency_opts['output'], 
                                    scaling=frequency_opts['scaling'],
                                    sides=frequency_opts['sides'], 
                                    taper=frequency_opts['taper'],
                                    halfbandwidth=frequency_opts['halfbandwidth'], 
                                    NW=frequency_opts['NW'],
                                    duration=frequency_opts['duration'], 
                                    overlap=frequency_opts['overlap'], 
                                    step=frequency_opts['step'],
                                    detrend=frequency_opts['detrend'], 
                                    n_fft=frequency_opts['n_fft'], 
                                    log_trans=frequency_opts['log_trans'], 
                                    centend=frequency_opts['centend'])
                    all_nSxx.append(Sxx)
                    
                    nSxx = ChanFreq()
                    nSxx.s_freq = Sxx.s_freq
                    nSxx.axis['freq'] = Sxx.axis['freq']
                    nSxx.axis['chan'] = Sxx.axis['chan']
                    nSxx.data = empty(1, dtype='O')
                    nSxx.data[0] = empty((Sxx.number_of('chan')[0],
                             Sxx.number_of('freq')[0]), dtype='f')
                    nSxx.data[0] = mean(
                            stack([x()[0] for x in all_nSxx], axis=2), axis=2)
            
            
            if filter_opts['laplacian']:
                segments.read_data(filter_opts['lapchan'], ref_chan)
            else:
                segments.read_data(chan, ref_chan)
            xfreq = []
            
            for sg, seg in enumerate(segments):
                print(f'Analysing segment {sg} of {len(segments)}')
                out = dict(seg)
                data = seg['data']
                timeline = data.axis['time'][0]
                out['start'] = timeline[0]
                out['end'] = timeline[-1]
                out['duration'] = len(timeline) / data.s_freq
                
                if frequency_opts['fast_fft']:
                    n_fft = next_fast_len(data.number_of('time')[0])
                else:
                    n_fft = frequency_opts['n_fft']
                
                if filter_opts['laplacian']:
                    selectchans = filter_opts['lapchan']
                else:
                    selectchans = chan
                
                if filter_opts['notch']:
                    print('Applying notch filtering.')
                    data.data[0] = notch_mne(data, oREF=filter_opts['oREF'], 
                                                channel=selectchans, 
                                                freq=filter_opts['notch_freq'],
                                                rename=filter_opts['laplacian_rename'],
                                                renames=filter_opts['renames'],
                                                montage=filter_opts['montage'])
                    
                if filter_opts['notch_harmonics']: 
                    print('Applying notch harmonics filtering.')
                    print(f'{selectchans}')
                    data.data[0] = notch_mne2(data, oREF=filter_opts['oREF'], 
                                              channel=selectchans, 
                                              rename=filter_opts['laplacian_rename'],
                                              renames=filter_opts['renames'],
                                              montage=filter_opts['montage'])    
                
                if filter_opts['bandpass']:
                    print('Applying bandpass filtering.')
                    data.data[0] = bandpass_mne(data, oREF=filter_opts['oREF'], 
                                              channel=selectchans,
                                              highpass=filter_opts['highpass'], 
                                              lowpass=filter_opts['lowpass'], 
                                              rename=filter_opts['laplacian_rename'],
                                              renames=filter_opts['renames'],
                                              montage=filter_opts['montage'])
                
                if filter_opts['laplacian']:
                    print('Applying Laplacian filtering.')
                    data.data[0] = laplacian_mne(data, 
                                         filter_opts['oREF'], 
                                         channel=chan, 
                                         ref_chan=ref_chan, 
                                         laplacian_rename=filter_opts['laplacian_rename'], 
                                         renames=filter_opts['renames'],
                                         montage=filter_opts['montage'])
                    data.axis['chan'][0] = asarray([x for x in chan])
                    selectchans = chan
                
                
                ### frequency transformation
                Sxx = frequency(data, output=frequency_opts['output'], 
                                scaling=frequency_opts['scaling'], 
                                sides=frequency_opts['sides'], 
                                taper=frequency_opts['taper'],
                                halfbandwidth=frequency_opts['halfbandwidth'], 
                                NW=frequency_opts['NW'],
                                duration=frequency_opts['duration'], 
                                overlap=frequency_opts['overlap'], 
                                step=frequency_opts['step'],
                                detrend=frequency_opts['detrend'], 
                                n_fft=n_fft, 
                                log_trans=frequency_opts['log_trans'], 
                                centend=frequency_opts['centend'])
                
                
                
                if norm:
        
                    for j, ch in enumerate(Sxx.axis['chan'][0]):
        
                        dat = Sxx.data[0][j,:]
                        sf = Sxx.axis['freq'][0]
                        f_res = sf[1] - sf[0] # frequency resolution
        
                        if norm == 'integral':
                            norm_dat = sum(dat) * f_res # integral by midpoint rule
                        else:
                            norm_dat = nSxx(chan=ch)[0]
        
                        Sxx.data[0][j,:] = dat / norm_dat
        
                out['data'] = Sxx
                
                
                if fooof_it:
                    
                    Fooofxx = frequency(data, output=frequency_opts['output'], 
                                    scaling=frequency_opts['scaling'], 
                                    sides=frequency_opts['sides'], 
                                    taper=frequency_opts['taper'],
                                    halfbandwidth=frequency_opts['halfbandwidth'], 
                                    NW=frequency_opts['NW'],
                                    duration=frequency_opts['duration'], 
                                    overlap=frequency_opts['overlap'], 
                                    step=frequency_opts['step'],
                                    detrend=frequency_opts['detrend'], 
                                    n_fft=n_fft, 
                                    log_trans=False, 
                                    centend=frequency_opts['centend'])
                    
                    freqs = Fooofxx.axis['freq'][0]
                    

                            
                    fooof_powers = zeros((len(chan), len(bands)))
                    fooof_ap_params = zeros((len(chan), len(bands), 2))
                    fooof_pk_params = ones((len(chan), len(bands), 9)) * nan
                    
                    for i in range(len(chan)):
                        fm.fit(freqs, Fooofxx.data[0][i], fooof_opts['freq_range'])
                        
                        for j, band in enumerate(bands):
                            fp = get_band_peak_fm(fm, band, fooof_opts['select_highest'],
                                                  threshold=fooof_opts['thresh_select'],
                                                  thresh_param=fooof_opts['thresh_param'])
                            if fp.ndim == 1:
                                fooof_powers[i, j] = gaussian_integral(fp[1], fp[2])
                            else:
                                pwr = asarray([gaussian_integral(fp[x, 1], fp[x, 2]) \
                                               for x in range(fp.shape[0])]).sum()
                                fooof_powers[i, j] = pwr
                                
                            # get fooof aperiodic parameters
                            fooof_ap_params[i, j, :] = fm.aperiodic_params_
                            
                            # get fooof peak parameters
                            fp = get_band_peak_fm(fm, band, False,
                                                  threshold=fooof_opts['thresh_select'],
                                                  thresh_param=fooof_opts['thresh_param'])
                            if fp.ndim == 1:
                                fooof_pk_params[i, j, :3] = fp
                            else:
                                n_peaks = min(fp.shape[0], 3)
                                fooof_pk_params[i, j, :n_peaks * 3] = fp[:n_peaks, 
                                                                      :].ravel()
                    
                    out['fooof_powers'] = fooof_powers
                    out['fooof_ap_params'] = fooof_ap_params
                    out['fooof_pk_params'] = fooof_pk_params
                            
                xfreq.append(out)
                
            if general_opts['freq_full']:                        
                if len(xfreq) == 1:
                    desc = None
                else:
                    as_matrix = asarray([y for x in xfreq for y in x['data']()[0]])
                    desc = get_descriptives(as_matrix)
                
                filename = outpath + p + '_' + vis + f'_freq_full{suffix}.csv'
                print('Writing to ' + filename)  
                export_freq(xfreq, filename, desc=desc)
                
            if general_opts['freq_band']:
                filename = outpath + p + '_' + vis + f'_freq_band{suffix}.csv'
                print('Writing to ' + filename)  
                export_freq_band(xfreq, frequency_opts['bands'], filename)
                
            if general_opts['freq_plot']:
                fig = Figure()
                FigureCanvas(fig)
                idx_plot = 1
                
                for seg in xfreq:
                    data = seg['data']
                    seg_chan = data.axis['chan'][0]
                    sf = data.axis['freq'][0]
                    if general_opts['max_freq_plot']:
                        idx_max = asarray(
                                [abs(x - general_opts['max_freq_plot']) for x in sf]).argmin()
                    for ch in seg_chan:
                        Sxx = data(chan=ch)[0]
                        ax = fig.add_subplot(len(xfreq), len(seg_chan), idx_plot)
                        ax.semilogy(sf[1:idx_max], Sxx[1:idx_max])
                        ax.set_xlabel('Frequency (Hz)')
                        
                        idx_plot += 1
                    
                print('Saving figure to ' + outpath + p + '.png')
                fig.savefig(outpath + p + suffix + '.png')
                
            if fooof_it:
                seg_info = ['Start time', 'End time', 'Duration', 'Stitches', 
                    'Stage', 'Cycle', 'Event type', 'Channel']
                #band_hdr = [str(b1) + '-' + str(b2) for b1, b2 in bands_fooof]
                band_hdr = [f'{b1}-{b2} Hz' for b1, b2 in bands]
                pk_params_hdr = ['peak1_CF', 'peak1_PW', 'peak1_BW', 
                              'peak2_CF', 'peak2_PW', 'peak2_BW', 
                              'peak3_CF', 'peak3_PW', 'peak3_BW', ]
                ap_params_hdr = ['Offset', 'Exponent']
                band_pk_params_hdr = ['_'.join((b, p)) for b in band_hdr 
                                      for p in pk_params_hdr]
                band_ap_params_hdr = ['_'.join((b, p)) for b in band_hdr 
                                      for p in ap_params_hdr]
                one_record = zeros((len(xfreq) * len(chan), 
                                    (len(seg_info) + len(bands) + 
                                    len(band_pk_params_hdr) + 
                                    len(band_ap_params_hdr))),
                                    dtype='O')
                for i, seg in enumerate(xfreq):
                    for j, ch in enumerate(chan):
                        cyc = None
                        if seg['cycle'] is not None:
                            cyc = seg['cycle'][2]
                        
                        one_record[i * len(chan) + j, :] = concatenate((asarray([
                            seg['start'],
                            seg['end'],
                            seg['duration'], 
                            seg['n_stitch'], # number of concatenated segments minus 1
                            seg['stage'],
                            cyc,
                            seg['name'], # event type
                            ch,
                            ]),
                            seg['fooof_powers'][j, :],
                            seg['fooof_pk_params'][j, ...].ravel(),
                            seg['fooof_ap_params'][j, ...].ravel(),
                            ))
                        
                outpath_fooof = outpath + p + '_' + vis + f'_fooofbands{suffix}.csv'
                print(f'Saving {outpath_fooof}')
                df = DataFrame(data=one_record, 
                               columns=(seg_info + band_hdr + band_pk_params_hdr
                                        + band_ap_params_hdr))
                df.to_csv(outpath_fooof)
    
    print("Job's done.")
    
    
    
def powerspec_summary_bands(root_dir, out_dir, part, visit, chan, stage, col_headers, 
                      excel_output_file):

    ## Script will create a file with multiple tab: once per STAGE*CHANNEL
    
    """ SCRIPT """
    csvs = []
    which_files = '*_freq_band.csv' # Name extension of output files to search for
    tabs = [' '.join((c,s)) for c in chan for s in stage] # Create matrices for each stage*channel combo
    
    # Loop through participants and visits to extract individual CSV output files
    if isinstance(part, list):
            None
    elif part == 'all':
            part = listdir(root_dir)
            part = [x for x in part if not('.' in x)]
    else:
        print("ERROR: 'part' must either be an array of subject ids or = 'all' ")      
            
    for i, p in enumerate(part):
            print(p)
            # loop through visits
            if visit == 'all':
                visit = listdir(root_dir + '/' + p)
                visit = [x for x in visit if not('.' in x)]
            
            for v, vis in enumerate(visit):
                ## Define files
                rdir = root_dir + p + '/' + vis + '/'
                
                try:
                    csvs.append(glob(f'{rdir}/{which_files}')[0]) #glob= va chercher tous les path du filename
                except Exception as e:
                    print(f'No files found for {p}, visit {vis}')
                
    ids = [x + '_' + y for x in part for y in visit] 
    idx_data_col = list(range(9,9+len(col_headers)))
    data = ones((len(chan) * len(stage), len(ids), len(idx_data_col))) * nan #create matrice - *nan to make sure that missing value will be blank
        
    for i, one_csv in enumerate(csvs): #on va passer dans la liste 1 a 1 (commence a compter a 0)    
        with open(one_csv, newline='') as f: #with open = to open a text/csv file
            csv_reader = reader(f, delimiter=',') #f = file, delimiter = what separate value in xcl (;, ou, ou .)
            
            for row in range(2): #skip first line
                row = next(csv_reader) #next = move forward
                
            if col_headers is None:
                col_headers = [row[x] for x in idx_data_col] #if None in col_header then will use the name of the csv
            
            for row in csv_reader:
                try: #try below, but if valueerror continue to except
                    int(row[0]) #find integer (chiffre)
                    idx_tab = tabs.index(' '.join((row[8], row[5]))) #join column 8 (channel) and 5(stage) to know which tab to go to
                    
                    for j, idx_col in enumerate(idx_data_col):
                        data[idx_tab, i, j] = row[idx_col] #create data, 1 column at a time and write it in correct tab
                                    
                except ValueError:
                    continue
    
    wb = Workbook()
    wb.save(excel_output_file)
                   
    with ExcelWriter(excel_output_file, engine="openpyxl", mode='a') as writer: #create xcl
    
        for i, tab_label in enumerate(tabs):
            df = DataFrame(data[i, ...], index=ids, columns=col_headers) #create dataframe and write data in it
            df.to_excel(writer, sheet_name=tab_label) #method of df to create xcl
        
    print(f'Saved to {excel_output_file}')
    
    
    
def powerspec_summary_full(root_dir, out_dir, part, visit, chan, stage, lowpass,
                      excel_output_file):

    ## Script will create a file with multiple tab: once per STAGE*CHANNEL
    
    """ SCRIPT """
    csvs = []
    which_files = '*_freq_full.csv' # Name extension of output files to search for
    tabs = [' '.join((c,s)) for c in chan for s in stage] # Create matrices for each stage*channel combo
    
    # Loop through participants and visits to extract individual CSV output files
    if isinstance(part, list):
            None
    elif part == 'all':
            part = listdir(root_dir)
            part = [x for x in part if not('.' in x)]
    else:
        print("ERROR: 'part' must either be an array of subject ids or = 'all' ")      
            
    for i, p in enumerate(part):
            print(p)
            # loop through visits
            if visit == 'all':
                visit = listdir(root_dir + '/' + p)
                visit = [x for x in visit if not('.' in x)]
            
            for v, vis in enumerate(visit):
                ## Define files
                rdir = root_dir + p + '/' + vis + '/'
                
                try:
                    csvs.append(glob(f'{rdir}/{which_files}')[0]) #glob= va chercher tous les path du filename
                except Exception as e:
                    print(f'No files found for {p}, visit {vis}')
                
    ids = [x + '_' + y for x in part for y in visit] 
    idx_data_col = list(range(9,9+lowpass*4))
    data = ones((len(chan) * len(stage), len(ids), len(idx_data_col))) * nan #create matrice - *nan to make sure that missing value will be blank
        
    for i, one_csv in enumerate(csvs): #on va passer dans la liste 1 a 1 (commence a compter a 0)    
        with open(one_csv, newline='') as f: #with open = to open a text/csv file
            csv_reader = reader(f, delimiter=',') #f = file, delimiter = what separate value in xcl (;, ou, ou .)
            
            for row in range(2): #skip first line
                row = next(csv_reader) #next = move forward
                
            
            col_headers = [row[x] for x in idx_data_col] #if None in col_header then will use the name of the csv
            
            for row in csv_reader:
                try: #try below, but if valueerror continue to except
                    int(row[0]) #find integer (chiffre)
                    idx_tab = tabs.index(' '.join((row[8], row[5]))) #join column 8 (channel) and 5(stage) to know which tab to go to
                    
                    for j, idx_col in enumerate(idx_data_col):
                        data[idx_tab, i, j] = row[idx_col] #create data, 1 column at a time and write it in correct tab
                                    
                except ValueError:
                    continue
    
    wb = Workbook()
    wb.save(excel_output_file)
                   
    with ExcelWriter(excel_output_file, engine="openpyxl", mode='a') as writer: #create xcl
    
        for i, tab_label in enumerate(tabs):
            df = DataFrame(data[i, ...], index=ids, columns=col_headers) #create dataframe and write data in it
            df.to_excel(writer, sheet_name=tab_label) #method of df to create xcl
        
    print(f'Saved to {excel_output_file}')