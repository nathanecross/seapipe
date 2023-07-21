#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 16:10:26 2021

@author: Nathan Cross
"""



from numpy import (asarray, float64, full, int64, nan)
from os import listdir, mkdir, path
from pandas import DataFrame, ExcelFile, read_csv
from wonambi.attr import Annotations


def sleepstats(in_dir,xml_dir,out_dir,rater,evt_type,stage,times,part='all',visit='all'):

    
    # 1. Sets lights off and on times for all recordings
    timesfile = ExcelFile(times).parse('Sheet1') #path in export_sleepstat
    
    if path.exists(out_dir):
                print(out_dir + " already exists")
    else:
        mkdir(out_dir)
    
    # 2. Set participants
    if isinstance(part, list):
        None
    elif part == 'all':
            part = listdir(in_dir)
            part = [ p for p in part if not '.' in p]
    else:
        print("ERROR: 'part' must either be an array of subject ids or = 'all' ")
       
    # 3. Run sleep statistic parameter extraction
    header = ['Visit', 'TIB_min', 'TotalWake_min', 'SL_min', 'WASOintra_min', 
              'Wmor_min', 'TSP_min', 'TST_min', 'SE', 'N1_min', 'N2_min', 
              'N3_min', 'REM_min', 'W_tsp', 'N1_tsp', 'N2_tsp', 'N3_tsp', 
              'REM_tsp', 'SSI', 'SFI', 'SL_toN2', 'SL_toN3', 'SL_toREM', 
              'SL_toNREM_5m', 'SL_toNREM_10m', 'SL_toN3_5m', 'SL_toN3_10m', 
              'Cycle_nbr', 'Arou_tot_h', 'Arou_N1_min', 'Arou_N2_min', 
              'Arou_N3_min', 'Arou_REM_min', 'Arou_NREM_min']
    nbrcol = len(header)-1
    
    if visit== 'all': #if multiple visits
            visit = listdir(in_dir + '/' + part[0])
            visit = [ x for x in visit if not '.' in x]
    
    sleepstatz = full((len(part), len(header)*len(visit)),nan, dtype=object) 

    sublist = [None] *len(part)
    part.sort()
    for i, p in enumerate(part): #for each participant...
        print(f'Running participant {p}')
        sublist[i] = p
        
        for v, vis in enumerate(visit): #for each visit...
            print(f'Running visit {vis}')   
            
            x_dir = xml_dir + p + '/' + vis + '/'
            
            if path.isdir(x_dir):
            
                xml_file = [x for x in listdir(x_dir) if x.endswith('.xml')]
            
                lights_off = timesfile[timesfile['ID'].str.contains(p)]
                lights_off = asarray(lights_off[lights_off['visit'].str.contains(vis)]['loff'])  
                if isinstance(lights_off[0],int64):
                    lights_off = float(lights_off[0])
                else:
                    lights_off = lights_off[0].astype(float64)
                lights_on = timesfile[timesfile['ID'].str.contains(p)]
                lights_on = asarray(lights_on[lights_on['visit'].str.contains(vis)]['lon'])
                if isinstance(lights_on[0],int64) or isinstance(lights_on[0],int):
                    lights_on = float(lights_on[0])
                else:
                    lights_on = lights_on[0].astype(float64)
                
                annot = Annotations(x_dir+xml_file[0], rater_name=rater)
                annot.export_sleep_stats(f'{out_dir}/{p}_{vis}_sleepstats.csv', lights_off, 
                                     lights_on)
                file = (out_dir + '/' + p + '_' + vis + '_sleepstats.csv')
                data = read_csv(file, sep=',', delimiter=None, header=1) #basic sleep stat of wonambi file
                ## create dataset sleep stat whole-night (26=nbr of column)
                sleepstatz[i,(nbrcol*v)+v+0] = vis
                sleepstatz[i,(nbrcol*v)+v+1] = data['Value 2'][4] # TIB (min)
                sleepstatz[i,(nbrcol*v)+v+2] = data['Value 2'][5] # Total Wake duration (min - between loff to lon)
                sleepstatz[i,(nbrcol*v)+v+3] = data['Value 2'][6] # Sleep latency (min - between loff to SO)
                sleepstatz[i,(nbrcol*v)+v+4] = data['Value 2'][7] # WASOintra (min - between SO to last epoch of sleep)
                sleepstatz[i,(nbrcol*v)+v+5] = data['Value 2'][8] # W morning (after last sleep epoch to Lon) (min)
                sleepstatz[i,(nbrcol*v)+v+6] = data['Value 2'][16] # TSP (min - between SO to last epoch of sleep)
                sleepstatz[i,(nbrcol*v)+v+7] = data['Value 2'][17] # TST (min - N1+N2+N3+REM)
                sleepstatz[i,(nbrcol*v)+v+8] = data['Value 1'][18] # Sleep efficiency (% - TST/TiB*100)
                sleepstatz[i,(nbrcol*v)+v+9] = data['Value 2'][10] # N1 (min)
                sleepstatz[i,(nbrcol*v)+v+10] = data['Value 2'][11] # N2 (min)
                sleepstatz[i,(nbrcol*v)+v+11] = data['Value 2'][12] # N3 (min)
                sleepstatz[i,(nbrcol*v)+v+12] = data['Value 2'][13] # REM (min)
                sleepstatz[i,(nbrcol*v)+v+13] = data['Value 1'][24] # Wake (% of TSP)  
                sleepstatz[i,(nbrcol*v)+v+14] = data['Value 1'][25] # N1 (% of TSP) # you can change  or add for %TST if necessary [29]
                sleepstatz[i,(nbrcol*v)+v+15] = data['Value 1'][26] # N2 (% of TSP) # you can change  or add for %TST if necessary [30]
                sleepstatz[i,(nbrcol*v)+v+16] = data['Value 1'][27] # N3 (% of TSP) # you can change  or add for %TST if necessary [31]
                sleepstatz[i,(nbrcol*v)+v+17] = data['Value 1'][28] # REM (% of TSP) # you can change  or add for %TST if necessary[32]
                sleepstatz[i,(nbrcol*v)+v+18] = data['Value 1'][35] # Stage shifts index (n/TSPhr)
                sleepstatz[i,(nbrcol*v)+v+19] = data['Value 1'][38] # Sleep fragmentation index (n/TSThr)
                sleepstatz[i,(nbrcol*v)+v+20] = data['Value 2'][40] # SL to N2 (min)
                sleepstatz[i,(nbrcol*v)+v+21] = data['Value 2'][41] # SL to N3 (min)
                sleepstatz[i,(nbrcol*v)+v+22] = data['Value 2'][42] # SL to REM (min)
                sleepstatz[i,(nbrcol*v)+v+23] = data['Value 2'][43] # SL to 5min of consecutive NREM (min)
                sleepstatz[i,(nbrcol*v)+v+24] = data['Value 2'][44] # SL to 10min of consecutive NREM (min)
                sleepstatz[i,(nbrcol*v)+v+25] = data['Value 2'][45] # SL to 5min of consecutive N3 (min)
                sleepstatz[i,(nbrcol*v)+v+26] = data['Value 2'][46] # SL to 10min of consecutive N3 (min)
                sleepstatz[i,(nbrcol*v)+v+27] = (len(annot.get_cycles())) # nbr of cycle  
                #sleepstatz[i,(nbrcol*v)+v+27] = 0 # nbr of cycle  #if no cycle
    
                ## add Arousal density
                file = (out_dir + '/' + p + '_' + vis +'_arou.csv')
                try:
                    annot.export_events(file, evt_type=evt_type,stage=stage) 
                except Exception as e:
                    annot.export_events(file, evt_type=evt_type) 
                data = read_csv(file, sep=',',skiprows=[0], delimiter=None)
                sleepstatz[i,(nbrcol*v)+v+28] = (len(data)/sleepstatz[i,7])*60 # Arousal index (n/hr)
                sleepstatz[i,(nbrcol*v)+v+29] = (len(data[data.Stage =='NREM1'])/sleepstatz[i,9]) #Arousal density N1 (n/min)
                sleepstatz[i,(nbrcol*v)+v+30] = (len(data[data.Stage =='NREM2'])/sleepstatz[i,10]) #Arousal density N2 (n/min)
                sleepstatz[i,(nbrcol*v)+v+31] = (len(data[data.Stage =='NREM3'])/sleepstatz[i,11]) #Arousal density N3 (n/min)
                sleepstatz[i,(nbrcol*v)+v+32] = (len(data[data.Stage =='REM'])/sleepstatz[i,12]) #Arousal density REM (n/min)  
                sleepstatz[i,(nbrcol*v)+v+33] = ((len(data[data.Stage =='NREM1'])+len(data[data.Stage =='NREM2'])+len(data[data.Stage =='NREM3']))
                                                 /(sleepstatz[i,9]+sleepstatz[i,10]+sleepstatz[i,11])) #Arousal density NREM (n/min) 
                
                
                if i==0 and v==0:
                    headerupdate = [h + '_' + vis for h in header]
                    #print(f'{headerupdate}')
                elif i==0:
                    headernew = [h + '_' + vis for h in header]
                    
                    headerupdate = headerupdate + headernew
                
            
            else:
                print(f'Participant {p} has no visit {vis}')
                if i==0 and v==0:
                    headerupdate = [h + '_' + vis for h in header]
                elif i==0:
                    headernew = [h + '_' + vis for h in header]
                    headerupdate = headerupdate + headernew
            
    
    # 3. Save sleep statistic parameters to file    

    keydf = DataFrame(sleepstatz, columns = headerupdate, index=sublist)
    DataFrame.to_csv(keydf, f'{out_dir}/sleepstats_all.csv', sep=',')


                  
    
    return    



def sleepstats_from_csvs(in_dir,out_dir,rater,stage,part='all',visit='all'):

    
    # 1. Sets lights off and on times for all recordings
    if path.exists(out_dir + '/sleepstats/'):
                print(out_dir + '/sleepstats/' + " already exists")
    else:
        mkdir(out_dir + '/sleepstats/')
    
    
    
    # 2. Set participants
    if isinstance(part, list):
        None
    elif part == 'all':
            files = [x for x in listdir(in_dir) if 'sleepstats.csv' in x]
            part = list({p.split('_')[0] for p in files})
    else:
        print("ERROR: 'part' must either be an array of subject ids or = 'all' ")
    
    part.sort()
        
    # 3. Run sleep statistic parameter extraction
    header = ['Visit', 'TIB_min', 'TotalWake_min', 'SL_min', 'WASOintra_min', 
              'Wmor_min', 'TSP_min', 'TST_min', 'SE', 'N1_min', 'N2_min', 
              'N3_min', 'REM_min', 'W_tsp', 'N1_tsp', 'N2_tsp', 'N3_tsp', 
              'REM_tsp', 'SSI', 'SFI', 'SL_toN2', 'SL_toN3', 'SL_toREM', 
              'SL_toNREM_5m', 'SL_toNREM_10m', 'SL_toN3_5m', 'SL_toN3_10m']
    nbrcol = len(header)-1
    
    if visit== 'all': #if multiple visits
            visit = list({p.split('_')[1] for p in files})
            
    visit.sort()
    
    sleepstatz = full((len(part), len(header)*len(visit)),nan, dtype=object) 

    sublist = [None]*len(part)
    end = 0
    for i, p in enumerate(part): #for each participant...
        print(f'Running participant {p}')
        print(f'i={i}')
        sublist[i] = p
        
        ivisit = [v.split('_')[1] for v in files if p in v]
        
        for v, vis in enumerate(ivisit): #for each visit...
            
            print(f'Running visit {vis}')   

            file = [f for f in files if p in f if vis in f]
            
            data = read_csv(in_dir + file[0], sep=',', delimiter=None, header=1) #basic sleep stat of wonambi file
            ## create dataset sleep stat whole-night (26=nbr of column)
            sleepstatz[i,(nbrcol*v)+v+0] = vis
            sleepstatz[i,(nbrcol*v)+v+1] = data['Value 2'][4] # TIB (min)
            sleepstatz[i,(nbrcol*v)+v+2] = data['Value 2'][5] # Total Wake duration (min - between loff to lon)
            sleepstatz[i,(nbrcol*v)+v+3] = data['Value 2'][6] # Sleep latency (min - between loff to SO)
            sleepstatz[i,(nbrcol*v)+v+4] = data['Value 2'][7] # WASOintra (min - between SO to last epoch of sleep)
            sleepstatz[i,(nbrcol*v)+v+5] = data['Value 2'][8] # W morning (after last sleep epoch to Lon) (min)
            sleepstatz[i,(nbrcol*v)+v+6] = data['Value 2'][16] # TSP (min - between SO to last epoch of sleep)
            sleepstatz[i,(nbrcol*v)+v+7] = data['Value 2'][17] # TST (min - N1+N2+N3+REM)
            sleepstatz[i,(nbrcol*v)+v+8] = data['Value 1'][18] # Sleep efficiency (% - TST/TiB*100)
            sleepstatz[i,(nbrcol*v)+v+9] = data['Value 2'][10] # N1 (min)
            sleepstatz[i,(nbrcol*v)+v+10] = data['Value 2'][11] # N2 (min)
            sleepstatz[i,(nbrcol*v)+v+11] = data['Value 2'][12] # N3 (min)
            sleepstatz[i,(nbrcol*v)+v+12] = data['Value 2'][13] # REM (min)
            sleepstatz[i,(nbrcol*v)+v+13] = data['Value 1'][24] # Wake (% of TSP)  
            sleepstatz[i,(nbrcol*v)+v+14] = data['Value 1'][25] # N1 (% of TSP) # you can change  or add for %TST if necessary [29]
            sleepstatz[i,(nbrcol*v)+v+15] = data['Value 1'][26] # N2 (% of TSP) # you can change  or add for %TST if necessary [30]
            sleepstatz[i,(nbrcol*v)+v+16] = data['Value 1'][27] # N3 (% of TSP) # you can change  or add for %TST if necessary [31]
            sleepstatz[i,(nbrcol*v)+v+17] = data['Value 1'][28] # REM (% of TSP) # you can change  or add for %TST if necessary[32]
            sleepstatz[i,(nbrcol*v)+v+18] = data['Value 1'][35] # Stage shifts index (n/TSPhr)
            sleepstatz[i,(nbrcol*v)+v+19] = data['Value 1'][38] # Sleep fragmentation index (n/TSThr)
            sleepstatz[i,(nbrcol*v)+v+20] = data['Value 2'][40] # SL to N2 (min)
            sleepstatz[i,(nbrcol*v)+v+21] = data['Value 2'][41] # SL to N3 (min)
            sleepstatz[i,(nbrcol*v)+v+22] = data['Value 2'][42] # SL to REM (min)
            sleepstatz[i,(nbrcol*v)+v+23] = data['Value 2'][43] # SL to 5min of consecutive NREM (min)
            sleepstatz[i,(nbrcol*v)+v+24] = data['Value 2'][44] # SL to 10min of consecutive NREM (min)
            sleepstatz[i,(nbrcol*v)+v+25] = data['Value 2'][45] # SL to 5min of consecutive N3 (min)
            sleepstatz[i,(nbrcol*v)+v+26] = data['Value 2'][46] # SL to 10min of consecutive N3 (min)
 
            #sleepstatz[i,(nbrcol*v)+v+27] = 0 # nbr of cycle  #if no cycle
            
            
            if i==0 and v==0:
                headerupdate = [h + '_' + vis for h in header]
                print(f'{headerupdate}')
            elif v>0 and end==0:
                end +=1
                headernew = [h + '_' + vis for h in header]
                headerupdate = headerupdate + headernew

    # 3. Save sleep statistic parameters to file    

    keydf = DataFrame(sleepstatz, columns = headerupdate, index=sublist)
    DataFrame.to_csv(keydf, f'{out_dir}/sleepstats_all.csv', sep=',')
             
    
    return   

