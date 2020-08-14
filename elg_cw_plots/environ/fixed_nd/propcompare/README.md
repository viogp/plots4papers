**mass_env_hist2.py** Compare stellar masses per environment.

Using box plots (https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.boxplot.html): 
line=median
box=25%,75%
whiskers=default is a bit random (1.5*quartiles values), set to [5,95] or use them to show the extent of the data (ax4.boxplot(data, showfliers=False)). 