import numpy as np 
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats.distributions import chi2
import warnings


samps = ['samp1A', 'samp1B', 'samp1C', 'samp2A', 'samp2B', 'samp2C']
#values = [0.600, 0.629, 0.728, 0.689, 0.523, 0.670]
values = [0.05, 0.1, 0.02, 0.9, 0.8, 0.82]

d = {}
d['Gene'] = ['thisgene'] * len(samps)
d['variable'] = samps
d['value'] = values
d['cond1'] = [1, 1, 1, 0, 0, 0]
d['cond2'] = [0, 0, 0, 1, 1, 1]
d['samples'] = [x + 1 for x in range(len(samps))]

rowdf = pd.DataFrame.from_dict(d)

with warnings.catch_warnings():
	warnings.filterwarnings('ignore')
	md = smf.mixedlm('value ~ cond1', data = rowdf, groups = 'cond1', missing = 'drop')
	mdf = md.fit(reml = False)
	nullmd = smf.mixedlm('value ~ 1', data = rowdf, groups = 'samples', missing = 'drop')
	nullmdf = nullmd.fit(reml = False)
	print mdf.llf
	print nullmdf.llf
	LR = 2 * (mdf.llf - nullmdf.llf)
	p = chi2.sf(LR, df = 1)
	print p
