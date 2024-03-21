#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:04:25 2022

@author: nathancross
"""

from os import listdir, mkdir, path
from . cfc_func import _allnight_ampbin
import matplotlib.pyplot as plt
from numpy import (append, array, argmax, arange, empty, floor, histogram, isnan,
                  mean, nanmean, nanstd, ndarray, pi, reshape, roll, sin, sqrt, std, zeros) 
from numpy.matlib import repmat
from safepickle import load
from pingouin import circ_mean
from scipy.stats import zscore




def plot_mean_amps(in_dir, out_dir, part, visit, chan, stage, band_pairs, 
                   colors, nbins, figargs={'dpi':400,'format':'svg'}):
    """
    Creates mean amplitude histograms across all events, per subject and visit, 
    to visualise the distribtion across the night for each subject.
    """
    
    plt.ioff()
    stagename = '_'.join(stage)
    
    # Check if band pairs is a list or a string
    if isinstance(band_pairs, str):
        band_pairs = [band_pairs]
    
    # Setup output directory
    if path.exists(out_dir):
            print(out_dir + " already exists")
    else:
        mkdir(out_dir)
    
    # Loop through band pairs    
    for i, bp in enumerate(band_pairs):
           
        # Set number of participants and visits
        if isinstance(part, list):
            None
        elif part == 'all':
                part = listdir(in_dir)
                part = [ p for p in part if not '.' in p]
        
        if isinstance(visit, list):
            None
        elif visit == 'all':
                visit = listdir(in_dir + part[0])
                visit = [ v for v in visit if not '.' in v]
        
        # Prepare output summaries
        datab = empty((len(part), len(visit), len(chan)),dtype='object')
        dist = zeros((len(part), len(visit), len(chan), 18))
        dist_sd = zeros((len(part), len(visit), len(chan), 18))
        w1 = zeros((len(part), len(visit), len(chan), 18))
        
        
        # Loop through channels
        for k, ch in enumerate(chan):
            
            # Prepare plots
            f, axarr = plt.subplots(nrows=len(part)+1, ncols=len(visit), 
                                    sharex=False, sharey=False, figsize=(18, 10),
                                    gridspec_kw={'height_ratios': append(repmat(3,1,len(part))[0],1)})
            f.suptitle(f'Mean_amps: {ch} - {bp}', fontsize=40)
            
            
            # Load in data
            #  Loop through participants
            for j,p in enumerate(part):
                
                # Loop through visits
                visit.sort()
                for v, vis in enumerate(visit):
                    
                    # Define files
                    rdir = in_dir + p + '/' + vis + '/'
                    ab_list = [x for x in listdir(rdir) if '.p' in x]
                    ab_file = [x for x in ab_list if bp in x if ch in x
                               if stagename in x]
                    ab_file = rdir + ab_file[0]
                
                    # Open files and extract data
                    with open(ab_file, 'rb') as file:
                        ab = load(file)
                        
                        dats = ab
                        datab[j,v,k] = dats
                        
                        dist[j,v,k,:] = nanmean(dats,axis=0) 
                        dist_sd[j,v,k,:] = nanstd(dats, axis=0) / sqrt(dats.shape[0])
                        w1[j,v,k,:], bin_edges = histogram(argmax(dats,axis=-1), 
                                                              bins=nbins) #range=(-pi, pi))
    
            
            # 1A.i. Mean amplitude histograms
            for x in range(0,len(part)):
                for y in range(0,len(visit)):
                    width = 2 * pi / len(dist[x,y,k,:])
                    pos = arange(width / 2, 2 * pi + width / 2, width)
                    
                    if not isinstance(axarr, ndarray): 
                        axarr.bar(pos, dist[x,y,k,:], width=width, color=colors[y], 
                                  edgecolor='w', yerr=dist_sd[x,y,k,:], ecolor='k', 
                                  capsize=0)
                        axarr.set_ylim([min(dist[x,y,k,:]-2*max(dist_sd[x,y,k,:])),
                                             max(dist[x,y,k,:])+2*max(dist_sd[x,y,k,:])])
                        axarr.yaxis.set_ticks([])
                        axarr.set_yticklabels([])
                        axarr.set_ylabel(part[x], fontsize=22)
                        if x == 0:
                                axarr.set_title(f'Visit: {visit[y]}', fontsize=32)
                    elif len(axarr.shape) <2:
                        axarr[x+y].bar(pos, dist[x,y,k,:], width=width, color=colors[y], 
                                       edgecolor='w', yerr=dist_sd[x,y,k,:], ecolor='k', 
                                       capsize=0)
                        axarr[x+y].set_ylim([min(dist[x,y,k,:]-2*max(dist_sd[x,y,k,:])),
                                             max(dist[x,y,k,:])+2*max(dist_sd[x,y,k,:])])
                        axarr[x+y].set_yticks([])
                        axarr[x+y].set_yticklabels([])
                        axarr[x+y].set_ylabel(part[x], fontsize=22)
                        if x == 0:
                                axarr[x+y].set_title(f'Visit: {visit[y]}', fontsize=32)
                    else:
                        axarr[x,y].bar(pos, dist[x,y,k,:], width=width, color=colors[y], 
                                       edgecolor='w', yerr=dist_sd[x,y,k,:], ecolor='k', 
                                       capsize=0)
                        axarr[x,y].set_ylim([min(dist[x,y,k,:]-2*max(dist_sd[x,y,k,:])),
                                             max(dist[x,y,k,:])+2*max(dist_sd[x,y,k,:])])
                        axarr[x,y].set_yticks([])
                        axarr[x,y].set_yticklabels([])
                        axarr[x,y].set_ylabel(part[x], fontsize=22)
                        if x == 0:
                                axarr[x,y].set_title(f'Visit: {visit[y]}', fontsize=32)
    
            # Last subplot (sine wave to illustrate phase)
            for sp in range(0,len(visit)):
                ax = axarr[-1,sp]
                x=arange(-pi, pi,pi/(nbins))
                ax.plot(x, sin(x))
                ax.set_xlim([-pi-0.2, pi+0.03])
                ax.set_xticks([-pi, -pi/2, 0, pi/2, pi])
                ax.set_xticklabels(['0', '', r'$\pi$', r'', r'2$\pi$'], fontsize=20)
                ax.set_ylim([-1.5,1.5])
                ax.set_yticklabels([])
                ax.set_yticks([])       
                    
            #print(f'bp = {type(bp)}; ch = {type(ch)}; stagename = {type(stagename)}')
            fname = bp + '_' + ch + '_' + stagename   
            f.savefig(f"{out_dir}/parts_hist_{fname}.{figargs['format']}", dpi=figargs['dpi'], 
                      format=figargs['format'])

        
def plot_prefphase(in_dir, out_dir, part, visit, chan, stage, band_pairs, 
                   colors, nbins, figargs={'dpi':400,'format':'svg'}):
    """
    Creates polar plots of preferred phase across all events, per subject and visit, 
    to visualise the distribtion across the night for each subject.
    """
    
    plt.ioff()
    stagename = '_'.join(stage)
    
    # Check if band pairs is a list or a string
    if isinstance(band_pairs, str):
        band_pairs = [band_pairs]
    
    # Setup output directory
    if path.exists(out_dir):
            print(out_dir + " already exists")
    else:
        mkdir(out_dir)
    
    # Loop through band pairs    
    for i, bp in enumerate(band_pairs):
           
        # Set number of participants and visits
        if isinstance(part, list):
            None
        elif part == 'all':
                part = listdir(in_dir)
                part = [ p for p in part if not '.' in p]
        
        if isinstance(visit, list):
            None
        elif visit == 'all':
                visit = listdir(in_dir + part[0])
                visit = [ v for v in visit if not '.' in v]
        
        # Prepare output summaries
        datab = empty((len(part), len(visit), len(chan)),dtype='object')
        dist = zeros((len(part), len(visit), len(chan), 18))
        dist_sd = zeros((len(part), len(visit), len(chan), 18))
        w1 = zeros((len(part), len(visit), len(chan), 18))
        
        
        # Loop through channels
        for k, ch in enumerate(chan):
            
            # Prepare plots
            f_p, axarr_p = plt.subplots(nrows=len(part), ncols=len(visit), 
                                        subplot_kw=dict(projection='polar'),
                                        figsize=(18, 10))
            f_p.suptitle(f'Preferred Phase: {ch} - {bp}', fontsize=40)
            
            
            # Load in data
            #  Loop through participants
            for j,p in enumerate(part):
                
                # Loop through visits
                visit.sort()
                for v, vis in enumerate(visit):
                    
                    # Define files
                    rdir = in_dir + p + '/' + vis + '/'
                    ab_list = [x for x in listdir(rdir) if '.p' in x]
                    ab_file = [x for x in ab_list if bp in x if ch in x
                               if stagename in x]
                    ab_file = rdir + ab_file[0]
                
                    # Open files and extract data
                    with open(ab_file, 'rb') as file:
                        ab = load(file)
                        
                        dats = ab
                        datab[j,v,k] = dats
                        
                        dist[j,v,k,:] = nanmean(dats,axis=0) 
                        dist_sd[j,v,k,:] = nanstd(dats, axis=0) / sqrt(dats.shape[0])
                        w1[j,v,k,:], bin_edges = histogram(argmax(dats,axis=-1), 
                                                              bins=nbins) #range=(-pi, pi))
    
                                
              # Preferred phase polar plots  
                for x in range(0,len(part)):
                    for y in range(0,len(visit)):
                        width = 2 * pi / len(dist[x,y,k,:])
                        if not isinstance(axarr_p, ndarray): 
                            axarr_p.bar(arange(0,2*pi,(2*pi)/18), w1[x,y,k,:], width=width, 
                                        color=colors[y], bottom=0.0)
                            axarr_p.xaxis.set_ticks([])
                            axarr_p.yaxis.set_ticks([])
                            axarr_p.set_xticklabels([])
                            axarr_p.set_yticklabels([])
                            axarr_p.set_ylabel(part[x], fontsize=22)
                            if x == 0:
                                    axarr_p.set_title(f'Visit: {visit[y]}', fontsize=32)
                        elif len(axarr_p.shape) <2:         
                            axarr_p[x+y].bar(arange(0,2*pi,(2*pi)/18), w1[x,y,k,:], 
                                             color=colors[y], width=width, bottom=0.0)
                            axarr_p[x+y].set_xticks([])
                            axarr_p[x+y].set_yticks([])
                            axarr_p[x+y].set_xticklabels([])
                            axarr_p[x+y].set_yticklabels([])
                            axarr_p[x+y].set_ylabel(part[x], fontsize=22)
                            if x == 0:
                                axarr_p[x+y].set_title(f'Visit: {visit[y]}', 
                                                            fontsize=32)
                        else:               
                            axarr_p[x,y].bar(arange(0,2*pi,(2*pi)/18), w1[x,y,k,:], 
                                             color=colors[y], width=width, bottom=0.0)
                            axarr_p[x,y].set_xticks([])
                            axarr_p[x,y].set_yticks([])
                            axarr_p[x,y].set_xticklabels([])
                            axarr_p[x,y].set_yticklabels([])
                            axarr_p[x,y].set_ylabel(part[x], fontsize=22)
                            if x == 0:
                                    axarr_p[x,y].set_title(f'Visit: {visit[y]}', 
                                                                fontsize=32)
           
            fname = bp + '_' + ch + '_' + stagename   
            f_p.savefig(f"{out_dir}/parts_polar_{fname}.{figargs['format']}", dpi=figargs['dpi'], 
                      format=figargs['format'])




def plot_prefphase_group(in_dir, out_dir, band_pairs, chan, cycle_idx, stage, nbins, 
                         norm, colors, comps = [('all','V1'), ('all','V2')], 
                         layout=None, figargs={'dpi':400,'format':'svg'}):
    
    '''
    This script generates a figure of the preferred phase from 
    
    comps = the comparisons you would like to plot, in the format
            of [(participants,visit)] 
            
            e.g. [('all', ['visit1']), 
                  ('all', ['visit2'])]
            
            or   [(['HC001','HC002'], ['visit1']),
                  (['PT001','PT002'], ['visit1'])]
            
    '''  
    # Setup output directory
    if path.exists(out_dir):
            print(out_dir + " already exists")
    else:
        mkdir(out_dir)
    
    # Check if band pairs is a list or a string
    if isinstance(band_pairs, str):
        band_pairs = [band_pairs]
        
    stagename = '_'.join(stage)
 
    # Loop through channels
    for k, ch in enumerate(chan):
        for b,bp in enumerate(band_pairs):
            print('')
            print("Plotting group-level preferred phase... ")
            print(f'CHANNEL: {ch}')
            print(f'BAND PAIR: {bp}')
            print('')
            
            # Prepare figures
            vecbin = zeros(nbins)
            width = 2 * pi / nbins
            for i in range(nbins):
                vecbin[i] = i * width + width / 2
            subplots = len(comps)
            if subplots<3:
                f_p, axarr_p = plt.subplots(1, 2, subplot_kw=dict(projection='polar'))
            elif subplots<5:
                f_p, axarr_p = plt.subplots(2, 2, subplot_kw=dict(projection='polar'))
            elif subplots<7:
                f_p, axarr_p = plt.subplots(3, 2, subplot_kw=dict(projection='polar'))
            else:
                if layout:
                    #f, axarr = plt.subplots(layout[0], layout[1], sharex=False, sharey='row')
                    f_p, axarr_p = plt.subplots(layout[0], layout[1], subplot_kw=dict(projection='polar'))
                else:
                    print("If plotting a large number (>6) of comparions, the parameter 'layout' must be specified.")
            
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
                for i, p in enumerate(part):
                    if visit == 'all':
                        visit = listdir(in_dir + '/' + p)
                        visit = [x for x in visit if not'.' in x]  
                    elif isinstance(visit,str):
                        visit = [visit]
                    
                # Define output object   
                datab = zeros((len(part),len(visit),nbins))
                
                # Loop through participants & visits
                for i, p in enumerate(part):
                    for j, vis in enumerate(visit): 
                        if not path.exists(in_dir + '/' + p + '/' + vis + '/'):
                            print(f'WARNING: input folder missing for Subject {p}, visit {vis}, skipping..')
                            continue
                        else:
                            p_files = [s for s in listdir(in_dir + '/' + p + '/' + vis) if bp in s if not s.startswith(".")]
                            p_files = [s for s in p_files if stagename in s] 
                            p_files = [s for s in p_files if ch in s]
                            
                            # Open files containing mean amplitudes (if existing)
                            if len(p_files) == 0:
                                print(f'WARNING: mean amplitudes file does not exist for Subject {p}, visit {j}, stage {stagename}, channel {ch} - check this. Skipping..')
                            elif len(p_files) >1:
                                print(f'WARNING: multiple mean amplitudes files exist for Subject {p}, visit {j}, stage {stagename}, channel {ch} - check this. Skipping..')
                            else:
                                print(f'Extracting... Subject {p}, visit {vis}')
                                ab_file = in_dir + '/' + p + '/' + vis + '/' + p_files[0]
                                with open(ab_file, 'rb') as f:
                                    ab = load(f)
                                
                                # Create bins for preferred phase
                                vecbin = zeros(nbins)
                                width = 2 * pi / nbins
                                for n in range(nbins):
                                    vecbin[n] = n * width + width / 2 
                                
                                
                                # Normalise mean amplitudes
                                try:
                                    ab = _allnight_ampbin(ab, 0, nbins, norm=norm)
                                except(ValueError):
                                    ab = ab / ab.sum(-1, keepdims=True)
                                    ab = ab.squeeze()
                                
                               

                                # Caculate mean direction for (binned) circular data
                                datab[i,j,:] = circ_mean(vecbin,mean(ab[:, :], axis=0) * 1000)
                                
                             
                w1, bin_edges = histogram(datab, bins=nbins, range=(-pi, pi))
                if subplots<3:
                    axarr_p[c].bar(bin_edges[:-1], w1, width=width, bottom=0.0, 
                                   color=colors[c])
                    axarr_p[c].set_xticks([])
                    axarr_p[c].set_xticklabels([])
                    axarr_p[c].title.set_text(comps[c])
                else:    
                    row = int(floor(subplots/c))
                    col = subplots % c
                    axarr_p[row,col].bar(bin_edges[:-1], w1, width=width, 
                                         bottom=0.0, color=colors[c])
                    axarr_p[row,col].set_xticks([])
                    axarr_p[row,col].set_xticklabels([])
                    axarr_p[row,col].title.set_text(comps[c])
                
            fname = bp + '_' + ch + '_' + stagename   
            f_p.savefig(f"{out_dir}/group_polar_{fname}.{figargs['format']}", dpi=figargs['dpi'], 
                      format=figargs['format'])
                
def plot_meanamps_group(in_dir, out_dir, band_pairs, chan, cycle_idx, stage, nbins, norm,
                         colors, comps = [('all','V1'), ('all','V2')], layout=None,
                        figargs={'dpi':400,'format':'svg'}):
    
    '''
    This script generates a figure of the group-level mean amps from 
    
    comps = the comparisons you would like to plot, in the format
            of [(participants,visit)] 
            
            e.g. [('all', ['visit1']), 
                  ('all', ['visit2'])]
            
            or   [(['HC001','HC002'], ['visit1']),
                  (['PT001','PT002'], ['visit1'])]
            
    '''  
    # Setup output directory
    if path.exists(out_dir):
            print(out_dir + " already exists")
    else:
        mkdir(out_dir)
    
    # Check if band pairs is a list or a string
    if isinstance(band_pairs,str):
        band_pairs = [band_pairs]
        
    stagename = '_'.join(stage)
 
    # Loop through channels
    for k, ch in enumerate(chan):
        for b,bp in enumerate(band_pairs):
            print('')
            print("Plotting group-level mean amps... ")
            print(f'CHANNEL: {ch}')
            print(f'BAND PAIR: {bp}')
            print('')
            
            # Prepare figures
            vecbin = zeros(nbins)
            width = 2 * pi / nbins
            for i in range(nbins):
                vecbin[i] = i * width + width / 2
            subplots = len(comps)
            if subplots<3:
                fig, axarr  = plt.subplots(2, 2, sharex=False, sharey=False,
                                        gridspec_kw={'height_ratios': [3, 1]})
            elif subplots<5:
                fig, axarr  = plt.subplots(3, 2, sharex=False, sharey=False,
                                        gridspec_kw={'height_ratios': [3, 3, 1]})
            elif subplots<7:
                fig, axarr  = plt.subplots(4, 2, sharex=False, sharey=False,
                                        gridspec_kw={'height_ratios': [3, 3, 3, 1]})
            else:
                if layout:
                    fig, axarr  = plt.subplots(layout[0]+1, layout[1], 
                                                sharex=False, sharey=False,
                                                gridspec_kw={'height_ratios': [repmat(3,1,layout[0]), 1]})
                else:
                    print("If plotting a large number (>6) of comparions, the parameter 'layout' must be specified.")
            
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
                for i, p in enumerate(part):
                    if visit == 'all':
                        visit = listdir(in_dir + '/' + p)
                        visit = [x for x in visit if not'.' in x]  
                    elif isinstance(visit,str):
                        visit = [visit]
                        
                # Define output object   
                datab = zeros((len(part),len(visit),nbins))
                
                # Loop through participants & visits
                for i, p in enumerate(part):
                    for j, vis in enumerate(visit): 
                        if not path.exists(in_dir + '/' + p + '/' + vis + '/'):
                            print(f'WARNING: input folder missing for Subject {p}, visit {vis}, skipping..')
                            continue
                        else:
                            p_files = [s for s in listdir(in_dir + '/' + p + '/' + vis) if bp in s if not s.startswith(".")]
                            p_files = [s for s in p_files if stagename in s] 
                            p_files = [s for s in p_files if ch in s]
                            
                            # Open files containing mean amplitudes (if existing)
                            if len(p_files) == 0:
                                print(f'WARNING: mean amplitudes file does not exist for Subject {p}, visit {j}, stage {stagename}, channel {ch} - check this. Skipping..')
                            elif len(p_files) >1:
                                print(f'WARNING: multiple mean amplitudes files exist for Subject {p}, visit {j}, stage {stagename}, channel {ch} - check this. Skipping..')
                            else:
                                print(f'Extracting... Subject {p}, visit {vis}')
                                ab_file = in_dir + '/' + p + '/' + vis + '/' + p_files[0]
                                with open(ab_file, 'rb') as f:
                                    ab = load(f)
                                
                                # Calculate mean amplitudes across night per subject
                                ab = nanmean(ab, axis=0)
                                
                               

                                # Caculate z-score for binned data
                                datab[i,j,:] = ab
                             
                                
                # Remove nans from output             
                databz = array([[[x if not isnan(x) else 0 for x in dim1] 
                                 for dim1 in dim2] for dim2 in datab])
                
                # Calculate mean and sd for comparison group
                data = reshape(databz,(databz.shape[0]*databz.shape[1],databz.shape[2]))
                data_m = nanmean(data,axis=0)
                data_sd = nanstd(data,axis=0) / sqrt(data.shape[0])
                
                # Prepare axes
                width = 2 * pi / nbins
                pos = arange(width / 2, 2 * pi + width / 2, width)
                
                if subplots<3:                
                    ax = axarr[0,c]
                else:
                    row = int(floor(subplots/c))
                    col = subplots % c
                    ax = axarr[row,col]
    
                # Set plotting limits
                lowlim = round(nanmean(data_m) - 2*nanmean(data_sd),2)
                sdlow = round(nanmean(data_m) - nanmean(data_sd),2)
                mid = round(nanmean(data_m),2)
                uplim = round(nanmean(data_m) + 2*nanmean(data_sd),2)
                sdup = round(nanmean(data_m) + nanmean(data_sd),2)
                
                # Plot data
                ax.bar(pos, data_m, width=width, color=colors[c], edgecolor='w', 
                            yerr=data_sd, ecolor='k', capsize=0)                
                ax.xaxis.set_ticks(pos)
                ax.set_xticklabels(arange(1,nbins+1,1).astype(str), fontsize='4')
                ax.yaxis.set_ticks([lowlim, sdlow, mid, sdup, uplim])
                ax.set_yticklabels([str(x) for x in 
                                    [lowlim, sdlow, mid, sdup, uplim]], fontsize='6')
                ax.set_ylim([lowlim, uplim])
                
                # Last subplot (sine wave to illustrate phase)
                for sp in range(0, 2):
                    ax = axarr[-1,sp]
                    x=arange(-pi, pi,pi/(nbins))
                    ax.plot(x, sin(x), color=colors[c])
                    ax.set_xlim([-pi-0.2, pi+0.03])
                    ax.set_xticks([-pi, -pi/2, 0, pi/2, pi])
                    ax.set_xticklabels(['0', '', r'$\pi$', r'', r'2$\pi$'], fontsize=10)
                    ax.set_ylim([-1.5,1.5])
                    ax.set_yticklabels([])
                    ax.set_yticks([])   
                
            fname = bp + '_' + ch + '_' + stagename   
            fig.savefig(f"{out_dir}/group_mean_amps_{fname}.{figargs['format']}", 
                        dpi=figargs['dpi'], 
                      format=figargs['format']) 
            
            
            