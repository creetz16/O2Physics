o2-analysis-lf-trackchecks --configuration json://config-tracks.json | 
o2-analysis-fdd-converter --configuration json://config-tracks.json |
o2-analysis-timestamp --configuration json://config-tracks.json | 
o2-analysis-trackextension --configuration json://config-tracks.json | 
o2-analysis-trackselection --configuration json://config-tracks.json | 
o2-analysis-event-selection --configuration json://config-tracks.json | 
o2-analysis-pid-tpc-full --add-qa 1 --configuration json://config-tracks.json | 
o2-analysis-multiplicity-table --configuration json://config-tracks.json > trackchecks_out.txt