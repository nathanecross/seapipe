
# Change Log

Beta Version

0.1 Initial release
0.2 All core functions included
0.3 Fixed dependencies when installing in fresh environments
0.3.1 Fixed a typo in dependencies
0.3.2 Added functionality of seabass (automatic sleep stage detection)
0.3.2 Bug fixes in stage and spindle detection for continuity when some subs crash
0.3.3 Fixed seabass when some channels (e.g. emg, eog) missing
0.3.4 More bug fixes in spindle detection for when running only some subjects at a time
0.3.5 Added functionality to read sub and ses from tracking sheet
0.3.7 Removed redundant scripts from pac
0.3.9 Added first artefact detection
0.4 Added continuous PAC and fixed bugs with event PAC
0.4.1 Updated format bids to allow subject-level exectution
0.4.2 Updates to export PAC summary
0.4.3 Improved loading channel names and completed PAC event dataset
0.4.4 Fixed bugs related to reading channel names from tracking sheet
0.4.5 Fixed bugs with power spectrum analysis and export
0.4.6 Fixed bugs with audit, whales, logs and added first REMS detection (YASA)
0.4.7 Corrected some more issues with logging, updates to artefact detection
0.4.8 Added artefact detection, fixed some loading issues, added back arousals to export macro stats
0.4.9 Updates to make_bids (added 'Woolcock' option and fixed bugs), whales
0.5.0 Updates to functions (make_bids, load, audit) and enabled data_dir as /sourcedata
0.5.1 Fixed bugs with channel renames (PSA and PAC) and make_bids
0.5.2 Fixed bugs, added function to infer polarity, and updated other (filtering) functions
0.5.3 Fixed bugs in PAC, misc and sleep staging