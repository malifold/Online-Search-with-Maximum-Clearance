In order to run the experiments shown for the star, run the file "starcomparisontime.py", execution takes a few seconds.


For experiments on the networks:
All relevant code is in "networks.py".
*testruns() will produce the graphs plotting clearance against time, for all 5 networks. Execution takes about an hour.

*berlinruns() will produce the 10 random runs we used. Execution takes about 10 minutes per run.

*chicagoruns() will produce the 45 random runs we used. Execution takes about one hour per run.


These runs for berlin and chicago will appear as pickled files in the folder "figs/cities". Once the runs for berlin or chicago have been produced, the functions
- compt(city)
- compr(city)
- compl(city)
plot results shown in Figures 10, 11, and 12 respectively, with city="ChicagoSketch". You can also use city="berlin" once the berlin runs are completed.

*competitiveness() will produce Figure 9.
