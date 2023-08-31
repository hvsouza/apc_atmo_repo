import uproot
import awkward as ak
import numpy as np
import pandas as pd
from tqdm import tqdm
import polars as pl
import os

class Anatree:
    tree:uproot.TTree
    fname: str

    def __init__(self, fname:str, entry_start=None, entry_stop=None, load_data = True):
        self.entry_start = entry_start
        self.entry_stop = entry_stop
        if fname.endswith('.root'):
            tree = uproot.open(f"{fname}:analysistree/anatree")
            self.tree = tree
        
        if load_data:
            self.load_anatree()

    def load_anatree(self, info=['nu', 'geant', 'reco_tracks', 'reco_showers']):
        if 'nu' in info: self._setup_nutree()
        if 'geant' in info: self._setup_geant()
        if 'reco_tracks' or 'reco_tracks' in info: self._setup_reco()

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

    def get_full_reco_tracks(self, df_tracks=None, df_geant=None, df_nu=None):
        if df_tracks is None:
            df_tracks = self.reco_tracks
        if df_geant is None:
            df_geant = self.geant
        if df_nu is None:
            df_nu = self.nu
        merged = df_tracks.join(df_geant, left_on=["subrun", "event", "trkg4id_pandoraTrack"], right_on=["subrun", "event", "TrackId_geant"], how="left")
        merged = merged.join(df_nu, left_on=["subrun", "event"], right_on=["subrun", "event"], how="inner")
        return merged
    
    def load_dqdx(self):
        subrun = self.tree['subrun'].array(library='np')
        event = self.tree['event'].array(library='np')
        ntracks = self.tree['ntracks_pandoraTrack'].array(library='np')
        dqdx = self.tree['trkdqdx_pandoraTrack'].array(library='ak')

        arr = {}

        arr['subrun'] = np.repeat(subrun, ntracks)
        arr['event'] = np.repeat(event, ntracks)

    def write_polars_parquet(self, fpath, batch_size = 10000, exclude=[''], max_events = None):
        """
        Function to create parquet files for polars. 
        fpath: str
            Path to save files
        batch_size: int
            Amount of events in each file created
        exclude: list
            Options: 'nu', 'geant' or 'reco'
        max_events : int
            Maximum number of events saved in total
        """ 

        # Setting up variables 
        from IPython.display import clear_output

        fpath = self._manage_path(fpath)

        exclude=list(exclude)
        files_to_write = ['nu', 'geant', 'reco']
        for exc in exclude:
            try:
                files_to_write.remove(exc)
            except:
                pass


        if 'reco' in files_to_write:
            files_to_write.remove('reco')
            files_to_write.append('reco_tracks')
            files_to_write.append('reco_showers')

        total_events = self.tree.num_entries
        if max_events != None:
            total_events = max_events
        nbatch = total_events//batch_size

        if total_events%batch_size!=0:
            nbatch+=1

        self.entry_start = 0
        self.entry_stop = batch_size # the `entry_stop` event is not executed 

        for file in files_to_write:
            file_to_test = f'{fpath}/{file}_0.parquet'
            file_exist = os.path.isfile(file_to_test)
            if file_exist:
                user_input = input(f'WARNING: There are parquet files that will be overwritten !\nContinue and overwrite it? (yes/no)?')
                yes_choices = ['yes', 'y']
                if user_input.lower() in yes_choices:
                    break
                else:
                    return

               
        output = []
        for b in range(nbatch):
            for out in output: print(out)
            print(f'Processing batch {b}')
            self.load_anatree()
            
            for file in files_to_write:
                if hasattr(self, file): # check if dataframe exist
                    df:pl.DataFrame
                    df = getattr(self, file)
                    print(f'Saving {file}...')
                    df.write_parquet(f'{fpath}/{file}_{b}.parquet')
                    
                else:
                    print(f"DataFrame '{file}' not found.")
            
            output.append(f'Batch {b+1}/{nbatch} done')
            clear_output()

            self.entry_start += batch_size
            self.entry_stop += batch_size

        for out in output: print(out)

    def _manage_path(self,fpath):
        if os.path.isfile(fpath):
            print('ERROR! `fpath` should be a path, not a file')
            assert(False)
        
        if os.path.exists(fpath):
            return fpath
        
        if os.path.isabs(fpath) == False:
            # Trying relative path
            
            relative_path = os.getcwd() + "/" + fpath
            if os.path.exists(relative_path):
                return relative_path
            else:
                fpath = relative_path
            
        user_input = input(f'Path {fpath} not found! Create one (yes/no)?')
        yes_choices = ['yes', 'y']
        # no_choices = ['no', 'n']
        if user_input.lower() in yes_choices:
            os.mkdir(fpath)
            return fpath
        else:
            assert(False)
    
    def read_parquet(self, fpath, batches=-1, batch_start = 0, types=['nu', 'geant', 'reco_tracks', 'reco_showers'], concat=True, **kwargs):
        '''
        Read parquet files
        fpath : str
            Path to files.
        batch : int
            Number of batches to be read. If -1, all will be read
        
        **kwargs
            Parameters for `read_parquet` of polars.

        
        types : str
            If you want to select what to load
        '''
        # self.nu = pl.DataFrame({})
        # self.geant = pl.DataFrame({})
        # self.reco_tracks = pl.DataFrame({})
        # self.reco_showers = pl.DataFrame({})
        if ~concat:
            self.nu = []
            self.geant = []
            self.reco_tracks = []
            self.reco_showers = []


        all_files = os.listdir(fpath)
        rechunk=False
        # loop by type
        for type in types:
            files = [f for f in all_files if type in f]
            files.sort()
            for b in range(batch_start):
                files.pop(0)
            if batches != -1:
                files = files[:-batches]

            nfiles = len(files)
            for i, file in enumerate(files):
                print(f'Reading {type} files... {i+1}/{nfiles}', end='\r')

                dfnew:pl.DataFrame
                dfnew = pl.scan_parquet(f'{fpath}/{file}', **kwargs)
                if hasattr(self, type): # check if dataframe exist
                    df:pl.DataFrame
                    df = getattr(self, type)
                    if concat:
                        if i == 0:
                            setattr(self,type,dfnew)
                        else:
                            if type == 'geant': rechunk = False
                            else: rechunk = True
                            df = pl.concat([df,dfnew], rechunk=rechunk)
                            setattr(self,type,df)
                    else:
                        df.append(dfnew)
                        setattr(self,type,df)
                else:
                    print(f"DataFrame '{file}' not found.")

                # if df_types[type].is_empty():
                # print('trying...')
                # print(self.df_types[type])
                # self.nu = pl.DataFrame({"4":2})
                # df_types[type] = pl.DataFrame({"2":1})

            print('')

                # else:
                #     df_types[type] = pl.concat([df_types[type],dfnew], rechunk=True)
            
                    
            
            
        


        


        


            

