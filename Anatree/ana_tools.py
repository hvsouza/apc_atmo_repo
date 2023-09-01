import matplotlib.pyplot as plt
import numpy as np
import polars as pl

def hstat(data, ax=None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    else:
        fig = plt.gcf()

    avg = np.nanmean(data.to_numpy())
    std = np.nanstd(data.to_numpy())

    if 'label' in kwargs:
        label = kwargs.pop('label')
    else:
        label = ""

    label = rf"{label} ($\mu={avg:.2f}$; $\sigma={std:.2f}$)"

    ax.hist(data, label=label, **kwargs)

    return fig

def get1Dcred(values, cred):
    #Assumes regular binning
    #Returns bins id that are within the specified credibility interval
    assert(0 < cred < 1)
    # bin_centers = 0.5*(bins[:-1] + bins[1:])
    bin_id = list(range(len(values)))
    zipped = list(zip(bin_id, values))
    bin_sorted = sorted(zipped, key= lambda x: x[1], reverse=True)

    id_sorted, values_sorted = zip(*bin_sorted)
    cum_ratio = np.cumsum(values_sorted)/np.sum(values_sorted)

    idx = np.searchsorted(cum_ratio, cred)
    return id_sorted[:idx]

def selection_events(extras = ['']):
    """
    Use this to quickly use subrun and event
    Ex: df.groupby(selection_events()).agg(
    ...
    )
    Or
    df.groupby(selection_events(['some_column','another'])).agg(
    ...
    )

    Also works with `pl.DataFrame.select`

    """
    r = ['subrun', 'event']
    if extras != ['']:
        r = r + extras
    return r

def get_event(subrun=0, event=1):
    """
    Quick function to get you a specific event
    """
    return (pl.col('subrun')==subrun) & (pl.col('event') == event)
