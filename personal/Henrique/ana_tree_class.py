from anytree import Node, RenderTree
import particle as pt
import uproot

class AnaTreeAnalyzer:
    def __init__(self, file = "", data = 0):
        self.file = file
        self.data = data
        self.tree = 0

        self.pdgs = 0
        self.ev = 0 
        self.nodes = 0

    def loadTree(self):
        self.tree = uproot.open(self.file)["anatree"]
        self.select_geant()
        self.select_truth()
        self.select_pandora()

        self.variables = self.geant_variables + self.truth_variables + self.pandora_variables
        self.data = self.tree.arrays(self.variables)

    def select_geant(self):
        self.geant_variables = [items for items in self.tree.keys() if 'geant' in items]
        self.geant_variables.remove('no_primaries_geant')
        self.geant_variables.remove('geant_list_size_geant')
        self.geant_variables.remove('geant_list_size_in_tpcAV_geant')
        self.geant_variables.remove('StartPointx_geant')
        self.geant_variables.remove('StartPointy_geant')
        self.geant_variables.remove('StartPointz_geant')
        self.geant_variables.remove('StartT_geant')
        self.geant_variables.remove('EndPointx_geant')
        self.geant_variables.remove('EndPointy_geant')
        self.geant_variables.remove('EndPointz_geant')
        self.geant_variables.remove('EndT_geant')
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

    def print_pdg_tree(self, mymaxlevel=2, **kwargs):
        parent_node = 0 # zero means that everyone will be printed
        nu_info = "enu_truth"
        try:
            parent_node = kwargs["node"]
            del kwargs["node"]
        except:
            parent_node = 0
        
        try:
            nu_info = kwargs["nu_info"]
            del kwargs["nu_info"]
        except:
            nu_info = "enu_truth"

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
            self.get_pdg_name()
            self.make_tree()
            self.print_pdg_tree(**kargs)
    

    # def get_node_by_pdg()



