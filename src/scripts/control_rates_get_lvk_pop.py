import paths
import json
import numpy as np
import pandas as pd

# Load LVK R&P posteriors for default PowerLaw+Peak model
with open(paths.rates_json, 'r') as f:
    rates_dict = json.load(f)

rates_df = pd.DataFrame(rates_dict['posterior']['content'])

# Rename entries because `sigma_chi` is actually a variance
rates_df.rename({'sigma_chi': 'var_chi'}, axis=1, inplace=True)
rates_df['log_prob'] = rates_df['log_likelihood'] + rates_df['log_prior']

# Identify MaP sample
i_maxp = np.argmax(rates_df['log_prob'].values)

# Compute Beta distribution alpha/beta parameters from its mean and variance
def get_alpha_beta(mu, var):
    # https://en.wikipedia.org/wiki/Beta_distribution#Mean_and_variance
    nu = mu*(1 - mu)/var - 1
    alpha = mu*nu
    beta = (1 - mu)*nu
    return alpha, beta

df = pd.DataFrame(np.array(get_alpha_beta(rates_df['mu_chi'],
                                          rates_df['var_chi'])).T,
                  columns=['alpha_chi', 'beta_chi'])
for k,v in df.items():
    rates_df[k] = v

rates_maxp = rates_df.iloc[i_maxp]
rates_maxp.to_csv(paths.rates_ref, sep=' ', header=False)
