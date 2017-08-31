These scripts and functions are used to to analyze spherical microindentatoin data using principles in Kalidindi and Pathak. Acta Materialia, 2008. and Pathak et al. Scipta Materialia. 2009 along with some of my own ideas.

The are many areas where these analysis can improved...

RunMe_Microindentation.m - use this to load, analyze, plot, etc. Most of the parameters that require editing are in this script. Use the raw .tra files for importing data.

MicroISS_Ph_v2.m - in here the stress-strain calculations are made including the zero point correction. For the zero point determination, see (hPnoCSM_v2.m). You may have to increase Prange and N for the load correction search.

hPnoCSM_v2.m - in here the load and dispalcement corrections are determined.

loadunload.m - this is used to determine the peaks and valleys of the load vs. time data.

No need to mess with these...
rsquare.m
mypolyfit.m

