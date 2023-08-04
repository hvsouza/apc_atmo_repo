# Anatree loader and tools

Here you can find the code to load the anatree root file in python. 
It uses [Uproot](https://uproot.readthedocs.io/en/latest/basic.html) and we use [Polars](https://pola-rs.github.io/polars-book/user-guide/) as the dataframe format. 

We also keep here some tools for plotting.

<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [Anatree loader and tools](#anatree-loader-and-tools)
- [How to use the loader](#how-to-use-the-loader)
- [Example of using polars and ana_tools](#example-of-using-polars-and-ana_tools)

<!-- markdown-toc end -->

# How to use the loader

Here is an example on how to use the class:  
Assuming you have this repository cloned on your `Documents` folder.  
You don't need to use `user_name`

``` python
# In jupyter-notebooks, this two comments allow to reload the class automatically in case of changes
%load_ext autoreload
%autoreload 2

import sys
import os

# using getlogin() returning username
user_name = os.getlogin()

sys.path.append(f'/home/{user_name}/Documents/Atmos_Pandora/apc_atmo_repo/Anatree/')
from anatree_class import Anatree
anatree:Anatree
```

And to load the data you can use the following

``` python
anatree = Anatree("/path/to/data/anatree_data.root, entry_stop=10000)
```

You have neutrino true information, Monte Carlo information, reconstructed tracks and showers stored  
in the following dataframe (df): `anatree.nu`, `anatree.geant`, `anatree.reco_tracks`, `anatree.reco_showers`

The df `anatree.merged` is the combination of reconstructed tracks with their respective MC particle information.

# Example of using polars and ana_tools

To plot the energy histogram of all primary muons, one can use the following:

``` python
import polars as pl

df = anatree.geant
df = df.filter(
    (pl.col('Mother_geant') == 0) & (pl.col('pdg_geant').abs() == 13)
    )

fig = hstat(df['Eng_geant'], label=r'Primary $\mu$', bins=200)
fig.legend()
ax = fig.gca()
ax.set_xlabel('Energy (GeV)')
ax.set_ylabel('Entries')
fig.show()
```



