import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from anytree import Node, RenderTree
nodes:Node

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


def merge_same_df(dataframes:list, filter_function=None, select_function=None):
    """
    Merged list of dataframe and return it as a dataframe

    Parameters
    ----------

    filter_function
        User defined function for filtering
    select_function
        User defined selection

    If nothing is passed, no filter is applying and all columns are selected

    """

    def filtering():
        if hasattr(filter_function, '__call__'):
            return filter_function()
        else:
            return True

    def selecting():
        if hasattr(select_function, '__call__'):
            return select_function()
        else:
            return pl.col('*')
        
    merged:pl.DataFrame
    merged = 0
    for i, df in enumerate(dataframes):
        df = df.filter(filtering()).collect()
        if i == 0:
            merged = df
        else:
            merged = pl.concat([merged,df])
    return merged


def make_tree(df:pl.DataFrame, subrun, event):
    nodes = {}
    def add_nodes(nodes, parent, child):
        if parent not in nodes:
            nodes[parent] = Node(parent)  
        if child not in nodes:
            nodes[child] = Node(child)
        nodes[child].parent = nodes[parent]
    df = df.filter(get_event(subrun=subrun,event=event))
    mother_and_selfid = df.select(['Mother_geant', 'TrackId_geant']).collect().to_numpy()
    for parent, child in zip(mother_and_selfid[:,0], mother_and_selfid[:,1]):
        add_nodes(nodes, parent, child)
        # print(nodes)
    return nodes


def print_tree(nu:pl.DataFrame, geant:pl.DataFrame, nodes:Node, subrun:int, event:int, maxlevel=2, particle_id = 0, only_ancestor=False, **kwargs):
    """
    Function to print nodes tree to see events. 
    
    Parameters
    ----------
    nu : pl.DataFrame
        Dataframe containing neutrino truth information

    geant : pl.DataFrame
        Dataframe containing geant information

    nodes : Node
        Nodes created with `make_tree`

    subrun, event : int
        Specific subrun and event

    maxlevel : int
        If printing entire tree, set the maxlevel to be printed.
        maxlevel=-1 will put no restriction

    particle_id : int
        If `particle_id != 0`, it will print the three starting for that specific particle

    only_ancestor : bool
        If set to `True`, only the parent of `particle_id` will be printed

    kwargs
        All the geant extra information desided. Example: `E="Eng_geant"`

    Examples
    --------
    As the `geant` dataframe from the anatree is usually very big, here is an example of how to use it

    >>> nodes = {}
    ... the_file = 0
    ... subruns = df.get_column('subrun').to_list()[:10]
    ... events = df.get_column('event').to_list()
    ... ids = df.get_column('trkg4id_pandoraTrack').to_list()
    ... for subrun, event, id in zip(subruns, events, ids):
    ...     for i, g in enumerate(anatree.geant):
    ...         nodes = make_tree(g, subrun, event)
    ...         if len(nodes) > 0:
    ...             the_file = i
    ...             break
    ...     print_tree( nu,
                        geant[the_file],
                        nodes,
                        particle_id=id,
                        subrun=subrun,
                        event=event,
                        only_ancestor=True,
                        maxlevel=3,
                        E='Eng_geant'
                        )
    0: nu(tau)~ E = 15.60 GeV 
    └── 3: K+ E = 4.13; 
        └── 3429: K+ E = 0.81; 
            └── 3494: mu+ E = 0.26; 

    """
    rendered = RenderTree(nodes[particle_id], maxlevel=maxlevel)
    if only_ancestor:
        try:
            nodes = nodes[particle_id].ancestors + (nodes[particle_id],) #add self as tuple 
            particle_id = 0
            rendered = []
            pre=''
            gap=''
            for i, node in enumerate(nodes):
                rendered.append([gap+pre,'none',node])
                gap=gap+'   '
                pre='└── '



        except IndexError:
            pass

    nu = nu.filter((get_event(subrun,event)))
    geant = geant.filter((get_event(subrun,event)))
    
    for pre, _, node in rendered:
        print(f"{pre}{node.name}:", end= " ")
        if node.name == 0:
            nupdg = particle.Particle.from_pdgid(nu.select('nuPDG_truth').item()).name
            print(fr"{nupdg} E = {nu.select('enu_truth').item():.2f} GeV", end=" ")
            print("")
        else:
            geant_info = geant.filter(pl.col('TrackId_geant')==node.name).select(
                pl.col('pdg_geant'),
                pl.col(list(kwargs.values()))
            ).collect()
            ppdg = geant_info.select('pdg_geant').item()
            # for key, val in kwargs.items():
                # print(f"{key}: {self.pdgs[node.name][1]}", end=" ")
            print(f"{particle.Particle.from_pdgid(ppdg)}", end=" ")
            for key, val in kwargs.items():
                extra_info = geant_info.select(pl.col(val)).item()
                print(f"{key} = {extra_info:.2f};", end=" ")
            print("")
