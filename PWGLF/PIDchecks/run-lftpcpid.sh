o2-analysis-lf-tpcpid --configuration json://config-lftpcpid.json | 
o2-analysis-timestamp --configuration json://config-lftpcpid.json | 
o2-analysis-track-propagation --configuration json://config-lftpcpid.json | 
o2-analysis-trackselection --configuration json://config-lftpcpid.json | 
o2-analysis-event-selection --configuration json://config-lftpcpid.json | 
o2-analysis-pid-tpc-full --add-qa 1 --configuration json://config-lftpcpid.json | 
o2-analysis-multiplicity-table --configuration json://config-lftpcpid.json > lftpcpid_out.txt