#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 16:10:26 2021

@author: Nathan Cross
"""

from datetime import datetime, date
from numpy import (asarray, float64, int64)
from os import listdir, path, walk
from pandas import DataFrame, read_csv
from wonambi.attr import Annotations
from ..utils.logs import create_logger, create_logger_outfile
                        
def export_sleepstats(xml_dir, out_dir, subs = 'all', sessions = 'all', 
               rater = None,  times = None, log_dir = None, 
               outfile = 'export_macro_stats_log.txt'):
    
    ### 0.a Set up logging
    flag = 0
    if not log_dir:
        log_dir = out_dir
    if outfile == True:
        today = date.today().strftime("%Y%m%d")
        now = datetime.now().strftime("%H:%M:%S")
        logfile = f'{log_dir}/export_sleep_macro_stats_{today}_log.txt'
        logger = create_logger_outfile(logfile=logfile, name='Export macro stats')
        logger.info('')
        logger.info(f"-------------- New call of 'Export macro stats' evoked at {now} --------------")
    elif outfile:
        logfile = f'{log_dir}/{outfile}'
        logger = create_logger_outfile(logfile=logfile, name='Export macro stats')
    else:
        logger = create_logger('Export macro stats')
    logger.info('')
    
    ## 1. Get lights on & lights off file
    if isinstance(times, str):
        times = read_csv(times, sep='\t')
    elif isinstance(times, DataFrame):
        None
    else:
        logger.critical('To export macro stats, lights off and lights on times are requried.')
        logger.info('Check documentation for how to export macro statistics:')
        logger.info('https://seapipe.readthedocs.io/en/latest/index.html')
        logger.info('-' * 10)
    
    # 2.a Get subjects
    if isinstance(subs, list):
        None
    elif subs == 'all':
        try:
            subs = next(walk(xml_dir))[1]
        except:
            logger.critical(f"{xml_dir} doesn't exist!") 
            return
        if len(subs) == 0:
            logger.critical(f'{xml_dir} is empty!') 
            return
    else:
        logger.info('')
        logger.critical("'subs' must either be an array of Participant IDs or = 'all' ")
        return
    subs.sort()
    
    # 2.b. Get sessions
    sub_ses = {}
    if isinstance(sessions, list):
        for sub in subs:
            sub_ses[sub] = sessions
    elif sessions == 'all':
        for sub in subs:
            session = next(walk(f'{xml_dir}/{sub}'))[1]
            session.sort()
            sub_ses[sub] = session
    else:
        logger.info('')
        logger.critical("'sessions' must either be an array of Session IDs or = 'all' ")
        return
    
    for s, sub in enumerate(sub_ses): #for each participant...
        if not sub_ses[sub]:
            logger.warning(f'No visits found in {xml_dir}/{sub}. Skipping..')
            flag +=1
            break
        for v, ses in enumerate(sub_ses[sub]): #for each visit...
            logger.debug(f'Export macro statistics for {sub}, {ses}')

            # 3.d. Get proper annotations file
            xml_file = [x for x in listdir(xml_dir +  '/' + sub + '/' + ses ) 
                        if x.endswith('.xml') if not x.startswith('.')]    
            if len(xml_file) == 0:                
                logger.warning(f'No sleep staging file found for {sub}, {ses}. Skipping..')
                flag+=1
                break
            elif len(xml_file) > 1:
                logger.warning(f"> 1 sleep staging file was found for {sub}, visit {ses} - to select the correct file you must define the variable 'keyword'. Skipping..")
                flag+=1
                break
            else:
                xml_file_path = f'{xml_dir}/{sub}/{ses}/{xml_file[0]}'
           
                # Search for sub and ses in tracking sheet
                l_times = times[times['sub']==sub]
                if l_times.size == 0:
                    logger.warning(f"Participant not found in column 'sub' in 'tracking.tsv' for {sub}, {ses}. Skipping...")
                    flag +=1
                    break
                l_times = l_times[l_times['ses']==ses]
                if l_times.size == 0:
                    logger.warning(f"Session not found in column 'ses' in 'tracking.tsv' for {sub}, {ses}. Skipping...")
                    flag +=1
                    break
                
                ## Lights OFF
                lights_off = l_times['loff']
                lights_off = asarray(lights_off.dropna())
                if lights_off.size == 0:
                    logger.warning(f"Lights Off time not found in 'tracking.tsv' for {sub}, {ses}. Skipping...")
                    flag +=1
                    break
                else:
                    if isinstance(lights_off[0],int64):
                        lights_off = float(lights_off[0])
                    else:
                        lights_off = lights_off.astype(float64)[0]
                
                ## Lights ON
                lights_on = l_times['lon']
                lights_on = asarray(lights_on.dropna())
                if lights_on.size == 0:
                    logger.warning(f"Lights On time not found in 'tracking.tsv' for {sub}, {ses}. Skipping...")
                    flag +=1
                    break
                else:
                    if isinstance(lights_on[0],int64):
                        lights_on = float(lights_on[0])
                    else:
                        lights_on = lights_on.astype(float64)[0]
                
                annot = Annotations(xml_file_path, rater_name=rater)
                annot.export_sleep_stats(f'{out_dir}/{sub}/{ses}/{sub}_{ses}_sleepstats.csv', 
                                         lights_off, lights_on)
                
                
    ### 3. Check completion status and print
    if flag == 0:
        logger.debug('Macro statistics export finished without ERROR.')  
    else:
        logger.warning('Macro statistics export finished with WARNINGS. See log for details.')
    
    return   


def sleepstats_from_csvs(xml_dir, out_dir, subs = 'all', sessions = 'all',
                         log_dir = None, outfile = 'macro_dataset_log.txt'):
    
    ### 0.a Set up logging
    flag = 0
    if not log_dir:
        log_dir = out_dir
    if outfile == True:
        today = date.today().strftime("%Y%m%d")
        now = datetime.now().strftime("%H:%M:%S")
        logfile = f'{log_dir}/export_sleep_macro_stats_{today}_log.txt'
        logger = create_logger_outfile(logfile=logfile, name='Export macro stats')
        logger.info('')
        logger.info(f"-------------- New call of 'Macro dataset' evoked at {now} --------------")
    elif outfile:
        logfile = f'{log_dir}/{outfile}'
        logger = create_logger_outfile(logfile=logfile, name='Export macro stats')
    else:
        logger = create_logger('Export macro stats')
    logger.info('')
    
    
    # 1. Get subjects
    if isinstance(subs, list):
        None
    elif subs == 'all':
            subs = next(walk(xml_dir))[1]
    else:
        logger.info('')
        logger.critical("'subs' must either be an array of Participant IDs or = 'all' ")
        return
    subs.sort()
    
    # 2. Get sessions
    if isinstance(sessions, list):
        None
    elif sessions == 'all':
        sessions = []
        for i, sub in enumerate(subs):
            sessions.append(next(walk(f'{xml_dir}/{sub}'))[1])
        sessions = sorted(set([x for y in sessions for x in y]))
    else:
        logger.info('')
        logger.critical("'sessions' must either be an array of Session IDs or = 'all' ")
        return
        
   
    
    # 3. Run sleep statistic parameter extraction
    header = ['ses', 'TIB_min', 'TotalWake_min', 'SL_min', 'WASOintra_min', 
              'Wmor_min', 'TSP_min', 'TST_min', 'SE_%', 'N1_min', 'N2_min', 
              'N3_min', 'REM_min', 'W_%tsp', 'N1_%tsp', 'N2_%tsp', 'N3_%tsp', 
              'REM_%tsp', 'SSI', 'SFI', 'SL_toN2_min', 'SL_toN3_min', 'SL_toREM_min', 
              'SL_toNREM_5m_min', 'SL_toNREM_10m_min', 'SL_toN3_5m_min', 'SL_toN3_10m_min']

    for v, ses in enumerate(sessions): #for each visit...
        df = DataFrame(data=None, index=subs, columns=header)
        
        for sub in subs: #for each participant...

            logger.debug(f'Extracting macro stats from {sub}, {ses}')  
            
            data_file = f'{xml_dir}/{sub}/{ses}/{sub}_{ses}_sleepstats.csv'
            
            if not path.exists(data_file):
                logger.warning(f"No (exported) staging data found for {sub}, {ses}. ")
                flag += 1
                continue
            else:
            
                data = read_csv(data_file, sep=',', delimiter=None, 
                            header=1, index_col=0) 
                df.loc['ses', sub] = ses
                df.loc['TIB_min', sub] = data['Value 2']['Total dark time (Time in bed)']
                df.loc['TotalWake_min', sub] = data['Value 2']['Wake duration']
                df.loc['SL_min', sub] = data['Value 2']['Sleep latency']
                df.loc['WASOintra_min', sub] = data['Value 2']['Wake after sleep onset']
                df.loc['Wmor_min', sub] = data['Value 2']['Wake, morning']
                df.loc['TSP_min', sub] = data['Value 2']['Total sleep period']
                df.loc['TST_min', sub] = data['Value 2']['Total sleep time']
                df.loc['SE_%', sub] = data['Value 1']['Sleep efficiency']
                df.loc['N1_min', sub] = data['Value 2']['N1 duration']
                df.loc['N2_min', sub] = data['Value 2']['N2 duration']
                df.loc['N3_min', sub] = data['Value 2']['N3 duration']
                df.loc['REM_min', sub] = data['Value 2']['REM duration']
                df.loc['W_%tsp', sub] = data['Value 1']['W % TSP']
                df.loc['N1_%tsp', sub] = data['Value 1']['N1 % TSP']
                df.loc['N2_%tsp', sub] = data['Value 1']['N2 % TSP']
                df.loc['N3_%tsp', sub] = data['Value 1']['N3 % TSP']
                df.loc['REM_%tsp', sub] = data['Value 1']['REM % TSP']
                df.loc['SSI', sub] = data['Value 1']['Switch index H']
                df.loc['SFI', sub] = data['Value 1']['Sleep fragmentation index H']
                df.loc['SL_toN2_min', sub] = data['Value 2']['Sleep latency to N2']
                df.loc['SL_toN3_min', sub] = data['Value 2']['Sleep latency to N3']
                df.loc['SL_toREM_min', sub] = data['Value 2']['Sleep latency to REM']
                df.loc['SL_toNREM_5m_min', sub] = data['Value 2']['Sleep latency to consolidated NREM, 5 min']
                df.loc['SL_toNREM_10m_min', sub] = data['Value 2']['Sleep latency to consolidated NREM, 10 min']
                df.loc['SL_toN3_5m_min', sub] = data['Value 2']['Sleep latency to consolidated N3, 5 min']
                df.loc['SL_toN3_10m_min', sub] = data['Value 2']['Sleep latency to consolidated N3, 10 min']

        df.to_csv(f'{out_dir}/{ses}_macro.csv', sep=',')
        
    ### 3. Check completion status and print
    logger.info('')
    if flag == 0:
        logger.debug('Macro statistics export finished without ERROR.')  
    else:
        logger.warning('Macro statistics export finished with WARNINGS. See log for details.')
    
    return  

