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

    def __init__(self, fname:str="", entry_start:int=None, entry_stop:int=None, load_data:bool = True, tree_name:str="analysistree/anatree"):
        self.entry_start = entry_start
        self.entry_stop = entry_stop
        if fname.endswith('.root'):
            tree = uproot.open(f"{fname}:{tree_name}")
            self.tree:uproot.TTree
            self.tree = tree
        
        if load_data:
            self.load_anatree()

    @classmethod
    def from_parquet(cls, *args, **kwargs) -> 'Anatree':
        obj = cls(load_data=False)
        obj.read_parquet(*args, **kwargs)
        return obj

    def load_anatree(self, info=['nu', 'geant', 'reco_tracks', 'reco_showers', 'pfp', 'reco_hits']):
        if 'nu' in info: self._setup_nutree()
        if 'geant' in info: self._setup_geant()
        if 'reco_tracks' in info: self._setup_reco_track()
        if 'reco_showers' in info: self._setup_reco_shower()
        if 'reco_hits' in info: self._setup_reco_hit()
        if 'pfp' in info: self._setup_pfp()

    def _setup_nutree(self):
        print("Loading nu infos")
        cols = ['run','subrun','event', 'nuPDG_truth','ccnc_truth', 'nuvtxx_truth','nuvtxy_truth','nuvtxz_truth',
                'enu_truth','nu_dcosx_truth','nu_dcosy_truth','nu_dcosz_truth',
                'lep_mom_truth', 'lep_dcosx_truth', 'lep_dcosy_truth', 'lep_dcosz_truth', 'mode_truth', 'nuWeight_truth',
                'Q2_truth', 'W_truth', 'X_truth', 'Y_truth', 'pot']
        
        ereco_cols = ['Ev_reco_nue', 'RecoLepEnNue', 'RecoHadEnNue', 'RecoMethodNue', 
                      'Ev_reco_numu', 'RecoLepEnNumu', 'RecoHadEnNumu', 'RecoMethodNumu', 
                      'LongestTrackContNumu', 'TrackMomMethodNumu', 
                      'RecoLepEnNumu_range', 'RecoLepEnNumu_mcs_chi2', 
                      'RecoLepEnNumu_mcs_llhd', 'Ev_reco_nc']
        cols = cols+ereco_cols
        arr = {}
        # nu_df = self.tree.arrays(cols, library='pd')
        # print(nu_df)
        for c in tqdm(cols):
            if c in self.tree.keys(filter_name=c):
                a = self.tree[c].array(library='ak', entry_start = self.entry_start, entry_stop = self.entry_stop)
                arr[c] = ak.ravel(a)

        nu = pd.DataFrame(arr)

        cols_reco = ['nuvtxx','nuvtxy','nuvtxz']
        if 'nnuvtx' in self.tree.keys(filter_name='nnuvtx'):
            nnuvtx = self.tree['nnuvtx'].array(library='ak', entry_start = self.entry_start, entry_stop = self.entry_stop)
        else:
            self.nu = pl.from_pandas(nu)
            return
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
        if 'geant_list_size_geant' in self.tree.keys(filter_name='geant_list_size_geant'):
            ngeant = self.tree['geant_list_size_geant'].array(library='ak', entry_start = self.entry_start, entry_stop = self.entry_stop)
        else:
            self.geant = pl.DataFrame()
            return
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

    def _setup_reco_shower(self):
        print("Loading shower infos")
        
        subrun = self.tree['subrun'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
        event = self.tree['event'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
   
        if not self.tree.keys(filter_name='nshowers_pandoraShower'):
            self.reco_showers = pl.DataFrame()
            return
        nshowers = self.tree['nshowers_pandoraShower'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)

        arr = {}

        arr['subrun'] = np.repeat(subrun, nshowers)
        arr['event'] = np.repeat(event, nshowers)

        cols = [key for key in self.tree.keys() if 'shwr' in key]
        cols.append('showerID_pandoraShower')
        for c in tqdm(cols):
            if c not in self.tree.keys(filter_name=c): continue
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
        if not self.tree.keys(filter_name='ntracks_pandoraTrack'):
            self.reco_tracks = pl.DataFrame()
            return
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
            if c not in self.tree.keys(filter_name=c): continue
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

    def _setup_reco_hit(self):
        print("Loading hit infos")
        run = self.tree['run'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
        subrun = self.tree['subrun'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
        event = self.tree['event'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
        
        nhits = self.tree['no_hits_stored'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)

        arr = {}

        arr['run'] = np.repeat(run, nhits)
        arr['subrun'] = np.repeat(subrun, nhits)
        arr['event'] = np.repeat(event, nhits)
        cols = self.tree.keys(filter_name="hit_*")

        for c in tqdm(cols):
            if c not in self.tree.keys(filter_name=c): continue
            a = self.tree[c].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
            flat = np.concatenate(a)
            arr[c] = flat
        self.reco_hits = pl.from_pandas(pd.DataFrame(arr))

    def _setup_pfp(self):
        print("Loading PFP infos")
        
        subrun = self.tree['subrun'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
        event = self.tree['event'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
   
        if 'nPFParticles' not in self.tree.keys(filter_name='nPFParticles'):
            self.pfp = pl.DataFrame()
            return
        npfps = self.tree['nPFParticles'].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)

        arr = {}

        arr['subrun'] = np.repeat(subrun, npfps)
        arr['event'] = np.repeat(event, npfps)

        cols = [key for key in self.tree.keys() if 'pfp' in key]
        cols.remove('pfp_daughterIDs')
        cols.remove('pfp_numClusters')
        cols.remove('pfp_clusterIDs')
        cols.remove('pfp_numNeutrinos')
        cols.remove('pfp_neutrinoIDs')
        for c in tqdm(cols):
            if c not in self.tree.keys(filter_name=c): continue
            a = self.tree[c].array(library='np', entry_start = self.entry_start, entry_stop = self.entry_stop)
            flat = np.concatenate(a)
            if 1 < len(flat.shape) <= 3:
                flat_x, flat_y, flat_z = np.split(flat, 3, axis=1)
                arr[f"{c}_x"] = flat_x.flatten()
                arr[f"{c}_y"] = flat_y.flatten()
                arr[f"{c}_z"] = flat_z.flatten()
            else:
                arr[c] = flat
            
        self.pfp = pl.from_pandas(pd.DataFrame(arr))

        temp_pfp = self.pfp
        temp_pfp = temp_pfp.join(self.pfp.groupby(['subrun','event']).agg(
            pl.col('pfp_selfID').filter( # Get ID of neutrino
                pl.col('pfp_isPrimary')==1
            ).first().alias('nuID')
        ),on=['subrun','event'],how='left').with_columns(
            pfp_parentID = pl.when( pl.col('pfp_parentID') == pl.col('nuID')).then(-1).otherwise(pl.col('pfp_parentID'))
        )
        
        temp = temp_pfp.select(
            pl.col('subrun'),
            pl.col('event'),
            pl.col('pfp_selfID'),
            pl.col('pfp_parentID'),
            pl.col('pfp_isTrack').alias('has_valid_pfp'),
        )
        if hasattr(self, 'reco_tracks') and not self.reco_tracks.is_empty():
            self.reco_tracks = self.reco_tracks.join(temp, left_on=['subrun','event','trkPFParticleID_pandoraTrack'], right_on=['subrun','event','pfp_selfID'], how='left')
        if hasattr(self, 'reco_showers') and not self.reco_showers.is_empty():
            self.reco_showers = self.reco_showers.join(temp, left_on=['subrun','event','shwr_PFParticleID_pandoraShower'], right_on=['subrun','event','pfp_selfID'], how='left')

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
    
    def get_full_pfp(self, df_pfp:pl.DataFrame=None, df_showers:pl.DataFrame=None, df_tracks:pl.DataFrame=None, df_geant:pl.DataFrame=None, df_nu:pl.DataFrame=None):
        if df_tracks is None:
            df_tracks = self.reco_tracks
        if df_showers is None:
            df_showers = self.reco_showers
        if df_pfp is None:
            df_pfp = self.reco_pfp
        if df_geant is None:
            df_geant = self.geant
        if df_nu is None:
            df_nu = self.nu

        merged = df_tracks.join(df_geant, left_on=["subrun", "event", "trkg4id_pandoraTrack"], right_on=["subrun", "event", "TrackId_geant"], how="left")
        merged = merged.join(df_nu, left_on=["subrun", "event"], right_on=["subrun", "event"], how="inner")

        merged = df_pfp.drop('pfp_parentID').join(
            merged,
            how='left',
            right_on=["subrun", "event", 'trkPFParticleID_pandoraTrack'],
            left_on=["subrun", "event", 'pfp_selfID']
            ).join(
            df_showers,
            how='left',
            right_on=["subrun", "event", 'shwr_PFParticleID_pandoraShower'],
            left_on=["subrun", "event", 'pfp_selfID']
        )

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
        files_to_write = ['nu', 'geant', 'reco', 'pfp']
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
            print(f'Processing batch {b+1}')
            self.load_anatree()
            
            for file in files_to_write:
                if hasattr(self, file): # check if dataframe exist
                    df:pl.DataFrame
                    df = getattr(self, file)
                    if df.is_empty():
                        print(f'Empty data {file} frame, not saving...')
                        continue
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
    
    def read_parquet(self, fpath, batches=-1, batch_start = 0, types=['nu', 'geant', 'reco_tracks', 'reco_showers', 'pfp'], concat=True, **kwargs):
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
        self.nu = []
        self.geant = []
        self.reco_tracks = []
        self.reco_showers = []
        self.pfp = []


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
                    print(f"DataFrame '{type}' not found.")

                # if df_types[type].is_empty():
                # print('trying...')
                # print(self.df_types[type])
                # self.nu = pl.DataFrame({"4":2})
                # df_types[type] = pl.DataFrame({"2":1})

            print('')

                # else:
                #     df_types[type] = pl.concat([df_types[type],dfnew], rechunk=True)
            
                    
            
            
        


        


        


            

