#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 13:55:32 2022

@author: nathancross
"""

from . cfc_func import mean_amp
from . mean_amps import generate_adap_bands
from numpy import (abs, argmax, asarray, ceil, concatenate, interp, linspace, 
                   logical_and, mean, median, nan, nanmean, newaxis, ones, pi, 
                   random, reshape, roll, sin, squeeze, vstack, where, zeros)
from os import listdir, mkdir, path, walk
from math import degrees, radians
from pandas import DataFrame
from pingouin import (circ_mean, circ_r, circ_rayleigh, circ_corrcc, circ_corrcl)
from shutil import copy
import sys
from tensorpac import Pac, EventRelatedPac
from termcolor import colored
from wonambi import Dataset
from wonambi.attr import Annotations
from wonambi.detect import match_events
from wonambi.trans import fetch




def pac_it_joint(rec_dir, xml_dir, out_dir, part, visit, cycle_idx, chan, rater, stage,
               polar, grp_name, cat, target, probe, buffer, ref_chan, oREF, nbins, idpac, 
               fpha, famp, dcomplex, filtcycle, width, min_dur, band_pairs, shift,
               laplacian=False, lapchan=None, adap_bands=(False,False),
               laplacian_rename=False, renames=None,
               filter_opts={'notch':False,'notch_harmonics':False, 'notch_freq':None,
                            'oREF':None,'chan_rename':False,'renames':None}):
   
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
    for p, pp in enumerate(part):
        if not path.exists(f'{out_dir}p/'):
            mkdir(f'{out_dir}p/')
        for v, vis in enumerate(visit):
            print('')
            print(f'Opening participant {pp}, visit {vis}')
            if not path.exists(f'{out_dir}p/{vis}/'):
                mkdir(f'{out_dir}p/{vis}/')
            out = zeros((len(part),10)) # intialise output dataframe
            for c, ch in enumerate(chan):
                print('')
                print(f'Analysing PAC for channel {ch}')
                chan_full = ch + ' ' + grp_name
                
                # Open edf and get information
                rec_file = [x for x in listdir(f'{rec_dir}/{pp}/{vis}/') if '.edf' in x]
                dataset = Dataset(f'{rec_dir}/{pp}/{vis}/{rec_file[0]}')
                s_freq = dataset.header['s_freq']
                nbuff = int(3*s_freq)
                
                # Open xml file and get information
                xml_file = [x for x in listdir(f'{xml_dir}/{pp}/{vis}/') if '.xml' in x]
                annot = Annotations(f'{xml_dir}/{pp}/{vis}/{xml_file[0]}')
                
                
                # EXTRACT BOTH EVENT TYPES
                # Events = (lower freq) phase signal
                events = [x for x in annot.get_events(target) if chan_full in x['chan']]
                events = events[1:]
                evt_mask = ones((len(events)), dtype=bool)
                
                # Probes = (higher freq) phase signal
                probes = [x for x in annot.get_events(probe) if chan_full in x['chan']]
                segs = fetch(dataset, annot, cat=(0,0,0,0), evt_type=[probe], cycle=None, 
                                 chan_full=[chan_full], reject_epoch=True, buffer=3,
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
                idx = [where(overlapping[x,:]==True)[0][0] for x in range(0,len(overlapping))]
                
                
                ## CALCULATE PREFERRED PHASE
                print('Extracting preferred phase of coupled events:')
                print(f'{target} : {probe}')
                print('')
                # Define frequency bands for PAC
                if adap_bands is True:
                    f_amp = famp[ch][pp + '_' + vis]
                    print(f'Using adapted bands for {p}, {vis}: {round(f_amp[0],2)}-{round(f_amp[1],2)} Hz')
                else:
                    f_amp = famp
                    print(f'Using fixed bands: freq={f_amp}')
                print('')
                
                # Define PAC object
                pac = Pac(idpac=idpac, f_pha=fpha, f_amp=f_amp, dcomplex='hilbert', 
                          cycle=(3,6), width=7, n_bins=nbins)
                
                # Prepare outputs
                pkbin = zeros((len(events),1))
                z=0 #(progress param)
                for e, ev in enumerate(events):
                    
                    # Print progress
                    z +=1
                    j = z/len(events)
                    sys.stdout.write('\r')
                    sys.stdout.write(f"Progress: [{'=' * int(50 * j):{50}s}] {int(100 * j)}%")
                    sys.stdout.write('\r')
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
                print('Extracting coupling strength.')
                print('')
                segs = fetch(dataset, annot, cat=(0,0,0,0), evt_type=[target], cycle=None, 
                                 chan_full=[chan_full], reject_epoch=True, buffer=3,
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
                    # Print progress
                    z +=1
                    j = z/len(segs_cp)
                    sys.stdout.write('\r')
                    sys.stdout.write(f"Progress: [{'=' * int(50 * j):{50}s}] {int(100 * j)}%")
                    sys.stdout.write('\r')
                    sys.stdout.flush()     
                    
                    # Extract data
                    dat = seg['data']()[0][0]
                    
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
                ab = roll(ampbin, shift, axis=-1)
                ma = nanmean(ab, axis=0)
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
            d = d.transpose()
            d.columns = ['mi_raw','mi_norm','mi_sig','pp_rad','pp_degrees','mvl',
                         'rho','pval','ray_z','ray_pv']
            stagename = ''.join(stage) 
            d.to_csv(path_or_buf=out_dir + '/' + p + '/' + vis + '/' + 
                     p + '_' + vis + '_' + ch + '_' + stagename + '_' + band_pairs + 
                     '_cfc_params.csv', sep=',')
    
    
    print('The function pac_it_joint completed without error.')  
              
                
                
                
                
                
                