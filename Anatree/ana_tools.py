import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import particle
import uproot
import pandas as pd
from tqdm import tqdm

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
    if not isinstance(extras, list):
        extras = [extras]
    r = ['subrun', 'event']
    if extras != ['']:
        r = r + extras
    print(r)
    return r

def get_event(subrun=0, event=1):
    """
    Quick function to get you a specific event
    """
    return (pl.col('subrun')==subrun) & (pl.col('event') == event)


def merge_same_df(dataframes:list, filter_function=None, select_function=None, df_concat = pl.DataFrame(), **kwargs):
    """
    Merged list of dataframe and return it as a dataframe

    Parameters
    ----------

    filter_function
        User defined function for filtering
    select_function
        User defined selection
    df_concat
        In case some selection is being done, df_concat is the dataframe to be concaternated
        Example of this is a dataframe with subrun and event to be used as filter
    kwargs
        Arguments to be passed to `join` with df_concat

    If nothing is passed for filter and selection, no filter is applying and all columns are selected
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
        df = df.filter(filtering()).select(selecting())
        if not df_concat.is_empty():
            df = df.join(df_concat.lazy(), **kwargs)
        df = df.collect()
        if i == 0:
            merged = df
        else:
            merged = pl.concat([merged,df])
    return merged

def check_fiducial(label:str, safe_x =20, safe_y=20, safe_z=20):
    limit_x = 363 - safe_x
    limit_y = 608 - safe_y
    upper_z = 1394 - safe_z
    lower_z = 0 + safe_z
    coord:str
    coord = ''
    coord_end = ''
    if label == 'pandoraTrack':
        coord = 'trkstart$_pandoraTrack'
        coord_end = coord.replace('start','end')
    elif 'pandoraShower' in label:
        coord = 'shwr_startx_pandoraShower'
        coord_end = coord.replace('start','end')
    elif 'geant' in label:
        coord = 'StartPoint$_'+label
        coord_end = 'EndPoint$_'+label
    else:
        return False
    return  (pl.col(coord.replace('$','x')).abs() < limit_x) & (pl.col(coord_end.replace('$','x')).abs() < limit_x) &\
            (pl.col(coord.replace('$','y')).abs() < limit_y) & (pl.col(coord_end.replace('$','y')).abs() < limit_y) &\
            (pl.col(coord.replace('$','z')).abs() < upper_z) & (pl.col(coord_end.replace('$','z')).abs() < upper_z) &\
            (pl.col(coord.replace('$','z')).abs() > lower_z) & (pl.col(coord_end.replace('$','z')).abs() > lower_z)



def make_tree(df:pl.LazyFrame, subrun, event):
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


def print_tree(nu:pl.DataFrame, geant:pl.LazyFrame, nodes:Node, subrun:int, event:int, maxlevel=2, particle_id = 0, only_ancestor=False, **kwargs):
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

def load_caf(file, save_it=False):
    rejections = ["Jagged", "Objects", "Group"]
    df:pl.DataFrame
    df = 0
    with uproot.open(file) as f:
        tree = f['cafTree']
        cols = [key for key in tree.keys() if not any([ x in str(tree[key].interpretation) for x in rejections])]
        data = f['cafTree'].arrays(cols, library='np')
    print('Converting to polars')
    df = pl.from_pandas(pd.DataFrame(data))
    for c in df.columns:
        df = df.rename({c: c.replace("/",".")})
    if save_it:
        file = file.replace('.root', '.parquet')
        df.write_parquet(file)
    return df

def proton_momentum_by_range(trkrange:pl.Expr) -> pl.Expr: 
    """
      Proton range-momentum tables from CSDA (Argon density = 1.4 g/cm^3):
      website: https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html

      CSDA values:
      double KE_MeV_P_Nist[31]={10, 15, 20, 30, 40, 80, 100, 150, 200, 250, 300,
      350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000,
      1500, 2000, 2500, 3000, 4000, 5000};

      double Range_gpercm_P_Nist[31]={1.887E-1,3.823E-1, 6.335E-1, 1.296,
      2.159, 7.375, 1.092E1, 2.215E1, 3.627E1, 5.282E1, 7.144E1,
      9.184E1, 1.138E2, 1.370E2, 1.614E2, 1.869E2, 2.132E2, 2.403E2,
      2.681E2, 2.965E2, 3.254E2, 3.548E2, 3.846E2, 4.148E2, 4.454E2,
      7.626E2, 1.090E3, 1.418E3, 1.745E3, 2.391E3, 3.022E3};

      Functions below are obtained by fitting power and polynomial fits to
      KE_MeV vs Range (cm) graph. A better fit was obtained by splitting the
      graph into two: Below range<=80cm,a a*(x^b) was a good fit; above 80cm, a
      polynomial of power 6 was a good fit

      Fit errors for future purposes:
      For power function fit: a=0.388873; and b=0.00347075
      Forpoly6 fit: p0 err=3.49729;p1 err=0.0487859; p2 err=0.000225834; p3
      err=4.45542E-7; p4 err=4.16428E-10; p5 err=1.81679E-13;p6
      err=2.96958E-17;

    ///////////////////////////////////////////////////////////////////////////
    //*********For proton, the calculations are valid up to 3.022E3 cm range
    //corresponding to a Muon KE of 5 GeV**********//
    ///////////////////////////////////////////////////////////////////////////
    """

    proton_mass = 938.272

    KE = pl.when(
        (trkrange > 0) & (trkrange < 80)
    ).then(
        29.9317*trkrange.pow(0.586304)
    ).otherwise(
        pl.when(
            (trkrange > 80) & (trkrange <= 3.022E3)
        ).then(
            149.904 + (3.34146 * trkrange) + (-0.00318856 * trkrange * trkrange) +
            (4.34587E-6 * trkrange * trkrange * trkrange) +
            (-3.18146E-9 * trkrange * trkrange * trkrange * trkrange) +
            (1.17854E-12 * trkrange * trkrange * trkrange * trkrange * trkrange) +
            (-1.71763E-16 * trkrange * trkrange * trkrange * trkrange * trkrange * trkrange)
        ).otherwise(
            0
        )
    )

    return (KE.pow(2) + 2*proton_mass*KE).sqrt()/1000

def manual_std(data:np.ndarray) -> float:
    a, b = np.quantile(data, [0.16, 0.84], method='linear')
    return (b - a)/2

def is_primary_reco() -> pl.Expr:
    return (pl.col('pfp_parentID') == -1) & (pl.col('pfp_isNeutrino') == 0)

def is_track() -> pl.Expr:
    return pl.col('pfp_isTrack') == 1
def is_shower() -> pl.Expr:
    return pl.col('pfp_isShower') == 1

def is_close_to_vertex(dist:float = 10) -> pl.Expr:
    return (
        (pl.col('trkstartx_pandoraTrack') - pl.col('nuvtxx'))**2 +
        (pl.col('trkstarty_pandoraTrack') - pl.col('nuvtxy'))**2 +
        (pl.col('trkstartz_pandoraTrack') - pl.col('nuvtxz'))**2
    ).sqrt() < dist

def trkke() -> pl.Expr:
    return pl.concat_list(
        pl.when(pl.col('trkke_pandoraTrack_x') > 0).then(pl.col('trkke_pandoraTrack_x')).otherwise(None),
        pl.when(pl.col('trkke_pandoraTrack_y') > 0).then(pl.col('trkke_pandoraTrack_y')).otherwise(None),
        pl.when(pl.col('trkke_pandoraTrack_z') > 0).then(pl.col('trkke_pandoraTrack_z')).otherwise(None)
    ).list.max()

def create_division(df:pl.DataFrame, particle, axes = ['x','y','z']):
    list_ndf = [f'trkpidndf_pandoraTrack_{ax}' for ax in axes]

    for ax in axes:
        name_of_div_col = f'div_{particle}_{ax}'
        df = df.select(
            pl.all(),
            (pl.col(f'trkpidchi{particle}_pandoraTrack_{ax}')/pl.col(f'trkpidndf_pandoraTrack_{ax}')).alias(name_of_div_col)
        ).with_columns(
            pl.when(pl.col(f'trkpidndf_pandoraTrack_{ax}') <0 ).then(None).otherwise(pl.col(name_of_div_col)).alias(name_of_div_col),
            pl.when(pl.col(f'trkpidndf_pandoraTrack_{ax}') <0 ).then(0).otherwise(1).alias(f'ndiv_{ax}'),
        )

    list_sum = [f'div_{particle}_{ax}' for ax in axes]
    list_div = [f'ndiv_{ax}' for ax in axes]

    df = df.with_columns( #first minimum, that works with nan
        (pl.min_horizontal(list_sum)).alias(f'trkpid{particle}_min'),
    ).fill_null(0).with_columns(
        (pl.sum_horizontal(list_sum)/pl.sum_horizontal(list_div)).alias(f'trkpid{particle}_av'),
        pl.max_horizontal(list_ndf).alias(f'trkpid{particle}_max')
       
    )
    return df

def pid_eval(df):
    df = create_division(df,'pr')
    df = create_division(df,'ka')
    df = create_division(df,'pi')
    df = create_division(df,'mu')
    return df

def particle_selection(df:pl.DataFrame, type='av'):
    """
    Make pid selection based on different methods
    type: list(str)
        `av` for average between chi2 x,y,z  \n
        `min` for getting minimum chi2/ndf  \n
        `max` for getting chi2 with maximum ndf  \n
    """

    df = df.with_columns(
        min_pid = pl.min(f'trkpidpr_{type}',f'trkpidka_{type}', f'trkpidpi_{type}',f'trkpidmu_{type}')
    ).with_columns(
        pid = pl.when(pl.col(f'trkpidpr_{type}')==pl.col('min_pid')).then(2212).otherwise(
        pl.when(pl.col(f'trkpidka_{type}')==pl.col('min_pid')).then(321).otherwise(
        pl.when(pl.col(f'trkpidpi_{type}')==pl.col('min_pid')).then(211).otherwise(
        pl.when(pl.col(f'trkpidmu_{type}')==pl.col('min_pid')).then(13).otherwise(
        0))))
    )

    return df

def manual_std(data:np.ndarray) -> float:
    a, b = np.quantile(data, [0.16, 0.84], method='linear')
    return (b - a)/2
