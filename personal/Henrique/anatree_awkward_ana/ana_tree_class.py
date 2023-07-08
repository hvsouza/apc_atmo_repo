from anytree import Node, RenderTree
import particle as pt
import uproot
import awkward as ak

class AnaTreeAnalyzer:
    tree:uproot.TTree
    nodes:Node
    def __init__(self, file = "", data = 0):
        self.file = file
        self.data = data
        self.tree = 0

        self.pdgs = 0
        self.ev = 0 
        self.nodes = 0
        self.pfp = []

    def loadTree(self):
        self.tree = uproot.open(self.file)
        self.select_general()
        self.select_geant()
        self.select_truth()
        self.select_pandora()
        self.select_pfp()

        self.create_alias()

        self.variables = self.general_variables + self.geant_variables + self.truth_variables + self.pandora_variables + self.pfp
        self.data = self.tree.arrays(self.variables+self.var_alias, aliases=self.alias)

    def create_alias(self):
        self.var_alias=["Length_geant","Length_pandoraTrack"]
        self.alias = { self.var_alias[0]: self._str_alias_geant(), self.var_alias[1]: self._str_alias_track() }

    def _str_alias_geant(self):
        return 'sqrt((EndPointx_drifted_geant - StartPointx_drifted_geant)**2'\
            '+ (EndPointy_drifted_geant - StartPointy_drifted_geant)**2' \
            '+ (EndPointz_drifted_geant - StartPointz_drifted_geant)**2)'

    def _str_alias_track(self):
        return 'sqrt((trkendx_pandoraTrack - trkstartx_pandoraTrack)**2 + '\
            '(trkendy_pandoraTrack - trkstarty_pandoraTrack)**2 + ' \
            '(trkendz_pandoraTrack - trkstartz_pandoraTrack)**2)'

    def select_general(self):
        self.general_variables = ['event', 'subrun']
        
    def select_pfp(self):
        self.pfp = [items for items in self.tree.keys() if 'pfp_' in items]
        self.pfp.remove('pfp_daughterIDs')
        self.pfp.remove('pfp_clusterIDs')
        self.pfp.remove('pfp_isNeutrino')
        self.pfp.remove('pfp_numNeutrinos')
        self.pfp.remove('pfp_neutrinoIDs')

    def select_geant(self):
        self.geant_variables = [items for items in self.tree.keys() if 'geant' in items]
        # self.geant_variables.remove('no_primaries_geant')
        # self.geant_variables.remove('geant_list_size_geant')
        # self.geant_variables.remove('geant_list_size_in_tpcAV_geant')
        self.geant_variables.remove('StartPointx_geant')
        self.geant_variables.remove('StartPointy_geant')
        self.geant_variables.remove('StartPointz_geant')
        # self.geant_variables.remove('StartT_geant')
        self.geant_variables.remove('EndPointx_geant')
        self.geant_variables.remove('EndPointy_geant')
        self.geant_variables.remove('EndPointz_geant')
        # self.geant_variables.remove('EndT_geant')
        self.geant_variables.remove("processname_geant")
        avtpc_variables = [items for items in self.tree.keys() if '_tpcAV_geant' in items]
        for val in avtpc_variables:
            # print(val)
            try:
                self.geant_variables.remove(val)
            except:
                pass

    def select_truth(self):
        self.truth_variables = [items for items in self.tree.keys() if 'truth' in items]

    def select_pandora(self):
        self.pandora_variables = [items for items in self.tree.keys() if 'pandora' in items]
        self.pandora_variables.remove('trkdedx_pandoraTrack')
        self.pandora_variables.remove('trkdqdx_pandoraTrack')
        self.pandora_variables.remove('trkresrg_pandoraTrack')
        self.pandora_variables.remove('trktpc_pandoraTrack')
        self.pandora_variables.remove('trkxyz_pandoraTrack')


    def create_pdg_ids(self):
        self.pdgs = {id: [ pdgcode, idx ] for idx,(id, pdgcode) in enumerate(zip(self.ev.TrackId_geant, self.ev.pdg_geant))}
        self.pdgs[0] = [self.ev["nuPDG_truth"][0], 0]

    def idx_to_id_geant(self, data = 0):
        try:
            if data==0:
                pass
        except:
                self.data = data
        idx_to_id = [{} for _,_ in enumerate(self.data)]

        for i, self.ev in enumerate(self.data):
            for idx, (id, pdgcode) in enumerate(zip(self.ev.TrackId_geant, self.ev.pdg_geant)):
                idx_to_id[i][id] = [pdgcode, idx]

        return idx_to_id

    def add_track_geant_info(self, data = 0):
        try:
            if data==0:
                pass
        except:
                self.data = data
        idx_to_id = self.idx_to_id_geant(self.data)

        tracks_to_g4_idx = []
        tracks_to_g4_pdg = []
        for ev, idx in zip(self.data, idx_to_id):
            if ev["ntracks_pandoraTrack"] != 0:
                pdg_id_of_track = ev["trkg4id_pandoraTrack"]
                pdg_idx = [idx[pdg_id][1] for pdg_id in pdg_id_of_track]
                pdg_pdg = [idx[pdg_id][0] for pdg_id in pdg_id_of_track]
                # print(ev["pdg_geant"][pdg_idx])
                # print(pdg_idx)
                tracks_to_g4_idx.append(pdg_idx)
                tracks_to_g4_pdg.append(pdg_pdg)
            else:
                # print("[]")
                tracks_to_g4_idx.append([])
                tracks_to_g4_pdg.append([])

        ak.from_regular(tracks_to_g4_idx)
        ak.from_regular(tracks_to_g4_pdg)

        self.data["trkg4idx_pandoraTrack"] = tracks_to_g4_idx
        self.data["trkg4pdg_pandoraTrack"] = tracks_to_g4_pdg

#data["pdg_id_to_idx"] = test


    def make_tree(self):
        def add_nodes(nodes, parent, child):
            if parent not in nodes:
                nodes[parent] = Node(parent)  
            if child not in nodes:
                nodes[child] = Node(child)
            nodes[child].parent = nodes[parent]

        self.nodes = {}  # store references to created nodes 
        for parent, child in zip(self.ev["Mother_geant"], self.ev["TrackId_geant"]):
            add_nodes(self.nodes, parent, child)
            # print(nodes)


    def get_pdg_name(self,pcode):
        try:
            ret = pt.Particle.from_pdgid(pcode).name
        except:
            ret = "nap"
        return ret

    def print_pdg_tree(self, mymaxlevel=2, nu_info = 'enu_truth', parent_node = 0, **kwargs):
        

        for pre, _, node in RenderTree(self.nodes[parent_node], maxlevel=mymaxlevel):
            print(f"{pre}{node.name}:{self.get_pdg_name(self.pdgs[node.name][0])}", end=" ")
            if node.name == 0:
                print(f"{self.ev[nu_info][0]:.4f}", end=" ")
                print("")
            else:
                for key, val in kwargs.items():
                    # print(f"{key}: {self.pdgs[node.name][1]}", end=" ")
                    print(f"{key}: {self.ev[val][self.pdgs[node.name][1]]:.4f}", end=" ")
                print("")

    def pdg_chain(self,range=[None,None], **kargs):
        for self.ev in self.data[range[0]:range[1]]:
            self.create_pdg_ids()
            self.make_tree()
            self.print_pdg_tree(3,**kargs)


    
    def check_fiducial(self, data, filter = None, type = "trk", xlim = 320, ylim = 560, zlow = 40, zhigh = 1350):
        if type == "trk":
            startx = "trkstartx_pandoraTrack"
            starty = "trkstarty_pandoraTrack"
            startz = "trkstartz_pandoraTrack"
            endx = "trkendx_pandoraTrack"
            endy = "trkendy_pandoraTrack"
            endz = "trkendz_pandoraTrack"

        elif type == "geant":
            startx = "StartPointx_drifted_geant"
            starty = "StartPointy_drifted_geant"
            startz = "StartPointz_drifted_geant"
            endx = "EndPointx_drifted_geant"
            endy = "EndPointy_drifted_geant"
            endz = "EndPointz_drifted_geant"

        try:
            if filter == None:
                filter = [True for _ in data]
        except:
            pass


        mask =  \
            (abs(data[startx][filter]) < xlim) &\
            (abs(data[endx][filter]) < xlim) &\
            (abs(data[starty][filter]) < ylim) &\
            (abs(data[endy][filter]) < ylim) &\
            (data[startz][filter] > zlow) &\
            (data[startz][filter] < zhigh) &\
            (data[endz][filter] > zlow) &\
            (data[endz][filter] < zhigh)

        return mask
    

    # def get_node_by_pdg()



