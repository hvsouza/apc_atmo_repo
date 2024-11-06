import uproot
import awkward as ak
import numpy as np
import pandas as pd
from tqdm import tqdm
import polars as pl
import os
import sys
user_name = os.getlogin()
sys.path.append(f'/home/{user_name}/Documents/Atmos_Pandora/apc_atmo_repo/Anatree/')
from  ana_tools import *



class PIDA:
    tree:uproot.TTree

    fname: str

    def __init__(self, fname:str="", entry_start:int=None, entry_stop:int=None, tree_name:str=":ana/t1"):

        self.entry_start = entry_start
        self.entry_stop = entry_stop
        self.fname = fname
        if not self.fname.endswith(".root"):
            self.fname += ".root"
        self.fname = self.fname + tree_name

        self.tree = uproot.open(self.fname)

        self.load_data()




    def get_singles(self, keys_single):
        print("Loading single values")
        arr = {}
        for c in tqdm(keys_single):
            a = self.tree[c].array(library='ak', entry_start = self.entry_start, entry_stop = self.entry_stop)
            arr[c] = ak.ravel(a)

        singles = pd.DataFrame(arr)

            
        self.dfsingles = pl.from_pandas(singles)
        print("Done...")

    def get_vectors(self, cols):
        print("Loading vectors")

        arr = {}

        arr['run'] = np.repeat(self.run, self.ntracks)
        arr['subrun'] = np.repeat(self.subrun, self.ntracks)
        arr['event'] = np.repeat(self.event, self.ntracks)

        for c in tqdm(cols):
            if c not in self.tree.keys(filter_name=c): continue
            if 'AsJagged' in str(self.tree[c].interpretation):
                a = self.tree[c].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
                flat = np.concatenate(a)
                arr[c] = flat
            else:
                a = self.tree[c].array(library='ak', entry_start = self.entry_start, entry_stop = self.entry_stop)
                aU = ak.Array([])
                aV = ak.Array([])
                aW = ak.Array([])
                aU = a[:,0,:]
                aV = a[:,1,:]
                aW = a[:,2,:]
                aU = ak.ravel(aU)
                aV = ak.ravel(aV)
                aW = ak.ravel(aW)

                arr[f"{c}_U"] = aU
                arr[f"{c}_V"] = aV
                arr[f"{c}_W"] = aW
        self.dfvector = pl.from_pandas(pd.DataFrame(arr))

        print("Done")

    def get_planes(self, cols):
        print("Loading planes")

        arr = {}

        arr['run'] = self.run
        arr['subrun'] = self.subrun
        arr['event'] = self.event

        for c in tqdm(cols):
            if c not in self.tree.keys(filter_name=c): continue
            a = self.tree[c].array(library='ak', entry_start = self.entry_start, entry_stop = self.entry_stop)
            aU = ak.Array([])
            aV = ak.Array([])
            aW = ak.Array([])
            aU = a[:,0]
            aV = a[:,1]
            aW = a[:,2]
            aU = ak.ravel(aU)
            aV = ak.ravel(aV)
            aW = ak.ravel(aW)

            arr[f"{c}_U"] = aU
            arr[f"{c}_V"] = aV
            arr[f"{c}_W"] = aW
        self.dfplanes = pl.from_pandas(pd.DataFrame(arr))

        print("Done")

    def load_data(self):

        allkeys = self.tree.keys()
        keys_single = []
        keys_vec = []
        keys_plane = ["allcalo_planes", "nallhits_planes"]
        exclusion = ['trknhitsmatch']

        for k in allkeys:
            interp = f"{self.tree[k].interpretation}"
            if k in keys_plane:
                pass
            elif ('AsObjects(AsVector' in interp or 'AsJagged' in interp):
                if k not in exclusion:
                    keys_vec.append(k)
            else:
                keys_single.append(k)
                
        self.run = self.tree['run'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
        self.subrun = self.tree['subrun'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
        self.event = self.tree['event'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
        # if 'ntracks' in allkeys:
        self.ntracks = self.tree['ntracks'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
        # else:
        #     if 'trkId' in allkeys:
        #         trkid = self.tree['trkId'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
        #         self.ntracks = [ len(tkid) for tkid in trkid]

        self.get_singles(keys_single)
        self.get_vectors(keys_vec)
        self.get_planes(keys_plane)
        # print(self.dfsingles)
        # print(self.dfvector)


def loaddata(path:str, entry_stop=None, forceit=False):
    parquet_path= path.replace(".root", ".parquet")
    if os.path.isfile(parquet_path) and not forceit:
        print(path, 'loaded parquet')
        df = pl.read_parquet(parquet_path)
        return df
    pid = PIDA(path, entry_stop=entry_stop)
    dfs:pl.DataFrame
    dfv:pl.DataFrame
    df:pl.DataFrame

    dfs = pid.dfsingles
    dfv = pid.dfvector
    dfp = pid.dfplanes

    df = dfs.join(dfv, on=selection_events())
    df = df.join(dfp, on=selection_events())
    df.write_parquet(file=parquet_path)
    return df