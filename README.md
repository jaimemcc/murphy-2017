# murphy-2017
Data files and analysis code for Murphy paper on protein preference

This repo provides access to all data and functions required to reproduce the figures in the manuscript, "Restriction of dietary protein leads to conditioned protein preference and elevated palatability of protein-containing food in rats" by Murphy et al submitted on 26 October 2017.

To run the code the following are required: 
(1) Med Associates data files
(2) Metafiles for behavioural data (cas9_metafile.txt) and body weight (cas9bw_metafile.txt)

These files are contained in a .zip file as supplementary material attached to the manuscript or can be downloaded from www.mccutcheonlab.com.

The code was prepared in Python 3.6 and most libraries should be included in a standard installation (e.g. via Anaconda). The exception is the statistics which require R and associated R-to-Python modules to be installed. To toggle statistics on and off change the variable statson which is assigned below. This variable is set to False by default so that code blocks that perform statistics will not run.

Any queries should be addressed to Dr James McCutcheon (jem64@le.ac.uk).
