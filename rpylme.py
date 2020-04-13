import gffutils
import os
import sys
import itertools
import gzip
import numpy as np
from Bio import SeqIO
import argparse
import subprocess
import pandas as pd
from collections import OrderedDict
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, Formula, FloatVector
from rpy2.rinterface import RRuntimeError
import math


pandas2ri.activate() #allow conversion between pandas dataframes and r dataframes

#define R packages
nlme = importr('nlme')
base = importr('base')
stats = importr('stats')
#qv = importr('qvalue')

#define formulae
fmla = Formula('value ~ 1 + cond1')
rndm = Formula('~ 1 | samples')
nullfmla = Formula('value ~ 1')
nullrndm = Formula('~1 | samples')

samps = ['samp1A', 'samp1B', 'samp1C', 'samp1D', 'samp2A', 'samp2B', 'samp2C', 'samp2D']
values = [0.1, 0.02, 0.05, 0.2, 0.8, 0.9, 0.85, 0.95]

d = {}
d['Gene'] = ['thisgene'] * len(samps)
d['variable'] = samps
d['value'] = values
d['cond1'] = [1, 1, 1, 1, 0, 0, 0, 0]
d['cond2'] = [0, 0, 0, 0, 1, 1, 1, 1]
d['samples'] = [x + 1 for x in range(len(samps))]

rowdf = pd.DataFrame.from_dict(d)

lm_alt = nlme.lme(fmla, random = rndm, data = rowdf, method = 'ML') #test
lm_null = nlme.lme(nullfmla, random = nullrndm, data = rowdf, method = 'ML') #control
logratio = (stats.logLik(lm_alt)[0] - stats.logLik(lm_null)[0]) * 2
pvalue = stats.pchisq(logratio, df = 1, lower_tail = False)[0]
#format decimal
pvalue = float('{:.2e}'.format(pvalue))

print rowdf
print pvalue