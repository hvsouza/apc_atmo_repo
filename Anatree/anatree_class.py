import uproot
import awkward as ak
import numpy as np
import pandas as pd
from tqdm import tqdm
import polars as pl

class Anatree:
    tree:uproot.TTree

    def __init__(self, fname, entry_start=None, entry_stop=None):
        self.entry_start = entry_start
        self.entry_stop = entry_stop
        tree = uproot.open(f"{fname}:analysistree/anatree")
        # tree.show()
        self.tree = tree


        self._setup_nutree()
        self._setup_geant()
        self._setup_reco()
        # self._process()

    def _setup_nutree(self):
        print("Loading nu infos")
        cols = ['nuPDG_truth','ccnc_truth','subrun','event','nuvtxx_truth','nuvtxy_truth','nuvtxz_truth',\
                'enu_truth','nu_dcosx_truth','nu_dcosy_truth','nu_dcosz_truth', 'mode_truth']
        arr = {}
        # nu_df = self.tree.arrays(cols, library='pd')
        # print(nu_df)
        for c in tqdm(cols):
            a = self.tree[c].array(library='ak', entry_start = self.entry_start, entry_stop = self.entry_stop)
            arr[c] = ak.ravel(a)

        nu = pd.DataFrame(arr)

        cols_reco = ['nuvtxx','nuvtxy','nuvtxz']
        nnuvtx = self.tree['nnuvtx'].array(library='ak', entry_start = self.entry_start, entry_stop = self.entry_stop)
        arr2 = {}
        arr2['subrun'] = np.repeat(arr['subrun'], nnuvtx)
        arr2['event'] = np.repeat(arr['event'], nnuvtx)
        
        for c in tqdm(cols_reco):
            a = self.tree[c].array(library='ak', entry_start = self.entry_start, entry_stop = self.entry_stop)
            arr2[c] = ak.ravel(a)

        nureco = pd.DataFrame(arr2)
        
        self.nu = pl.from_pandas(pd.merge(nu, nureco, how='left', on=['subrun', 'event']))
        
    def _setup_geant(self):
        print("Loading geant infos")
        subrun = self.tree['subrun'].array(library='ak', entry_start = self.entry_start, entry_stop = self.entry_stop)
        event = self.tree['event'].array(library='ak', entry_start = self.entry_start, entry_stop = self.entry_stop)
        ngeant = self.tree['geant_list_size_geant'].array(library='ak', entry_start = self.entry_start, entry_stop = self.entry_stop)

        arr = {}

        arr['subrun'] = np.repeat(subrun, ngeant)
        arr['event'] = np.repeat(event, ngeant)

        # idx = pd.MultiIndex.from_arrays([[subrun[i] for i in range(len(ngeant)) for _ in range(ngeant[i])], [event[i] for i in range(len(ngeant)) for j in range(ngeant[i])]], names=("subrun", "ebent"))

        cols = [key for key in self.tree.keys() if 'geant' in key]
        cols.remove("no_primaries_geant")
        cols.remove("geant_list_size_geant")
        cols.remove("geant_list_size_in_tpcAV_geant")
        
        for c in tqdm(cols):
            a = self.tree[c].array(library='ak', entry_start = self.entry_start, entry_stop = self.entry_stop)
            # arr[c] = np.concatenate(a)
            arr[c] = ak.ravel(a)

        self.geant = pl.from_pandas(pd.DataFrame(arr))

    def _setup_reco(self):
        self._setup_reco_shower()
        self._setup_reco_track()

    def _setup_reco_shower(self):
        print("Loading shower infos")
        
        subrun = self.tree['subrun'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
        event = self.tree['event'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
   
        nshowers = self.tree['nshowers_pandoraShower'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)

        arr = {}

        arr['subrun'] = np.repeat(subrun, nshowers)
        arr['event'] = np.repeat(event, nshowers)

        cols = [key for key in self.tree.keys() if 'shwr' in key]
        cols.append('showerID_pandoraShower')
        for c in tqdm(cols):
            a = self.tree[c].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
            flat = np.concatenate(a)
            if len(flat.shape) > 1:
                flat_x, flat_y, flat_z = np.split(flat, 3, axis=1)
                arr[f"{c}_x"] = flat_x.flatten()
                arr[f"{c}_y"] = flat_y.flatten()
                arr[f"{c}_z"] = flat_z.flatten()
            else:
                arr[c] = flat
            
        self.reco_showers = pl.from_pandas(pd.DataFrame(arr))

    def _setup_reco_track(self):
        print("Loading track infos")
        subrun = self.tree['subrun'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
        event = self.tree['event'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
        ntracks = self.tree['ntracks_pandoraTrack'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)

        arr = {}

        arr['subrun'] = np.repeat(subrun, ntracks)
        arr['event'] = np.repeat(event, ntracks)
        cols = [key for key in self.tree.keys() if 'trk' in key]

        exclusion_list = ['cosmic', 'T0', 'momms', 'pidmva', 'evtxid']  #all this are empty
        cols = [x for x in cols if all(exclusion_item not in x for exclusion_item in exclusion_list)]
        cols.append('trkId_pandoraTrack')
        cols.remove("trkdedx_pandoraTrack")
        cols.remove("trkdqdx_pandoraTrack")
        cols.remove("trkresrg_pandoraTrack")
        cols.remove("trktpc_pandoraTrack")
        cols.remove("trkxyz_pandoraTrack")

        for c in tqdm(cols):
            a = self.tree[c].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
            flat = np.concatenate(a)
            if len(flat.shape) > 1:
                flat_x, flat_y, flat_z = np.split(flat, 3, axis=1)
                arr[f"{c}_x"] = flat_x.flatten()
                arr[f"{c}_y"] = flat_y.flatten()
                arr[f"{c}_z"] = flat_z.flatten()
            else:
                arr[c] = flat
        self.reco_tracks = pl.from_pandas(pd.DataFrame(arr))

    def get_full_reco_tracks(self):
        merged = self.reco_tracks.join(self.geant, left_on=["subrun", "event", "trkg4id_pandoraTrack"], right_on=["subrun", "event", "TrackId_geant"], how="left")
        merged = merged.join(self.nu, left_on=["subrun", "event"], right_on=["subrun", "event"], how="inner")
        return merged
    
    def load_dqdx(self):
        subrun = self.tree['subrun'].array(library='np')
        event = self.tree['event'].array(library='np')
        ntracks = self.tree['ntracks_pandoraTrack'].array(library='np')
        dqdx = self.tree['trkdqdx_pandoraTrack'].array(library='ak')

        arr = {}

        arr['subrun'] = np.repeat(subrun, ntracks)
        arr['event'] = np.repeat(event, ntracks)

