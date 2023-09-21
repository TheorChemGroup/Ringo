import networkx as nx
from networkx.algorithms import isomorphism
import numpy as np
from numpy.linalg import norm

RADII = {'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.69, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76, 'Sc': 1.7, 'Ti': 1.6, 'V': 1.53, 'Cr': 1.39, 'Mn': 1.61, 'Fe': 1.52, 'Co': 1.5, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 'Ga': 1.22, 'Ge': 1.2, 'As': 1.19, 'Se': 1.2, 'Br': 1.2, 'Kr': 1.16, 'Rb': 2.2, 'Sr': 1.95, 'Y': 1.9, 'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.4, 'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89, 'Tm': 1.9, 'Yb': 1.87, 'Lu': 1.87, 'Hf': 1.75, 'Ta': 1.7, 'W': 1.62, 'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.4, 'At': 1.5, 'Rn': 1.5, 'Fr': 2.6, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06, 'Pa': 2.0, 'U': 1.96, 'Np': 1.9, 'Pu': 1.87, 'Am': 1.8, 'Cm': 1.69}


def write_sdf(G, p, index, fname):
    structure = p[index].xyz
    atom_symbols = p.atom_symbols
    lines = ["", "", ""]
    lines.append("%3d%3d  0  0  0  0  0  0  0  0999 V2000" % (G.number_of_nodes(), G.number_of_edges()))
    #   -5.5250    1.6470    1.0014 C   0  0  0  0  0  0  0  0  0  0  0  0
    new_indexing = list(G.nodes)
    for atom in new_indexing:
        lines.append("%10.4f%10.4f%10.4f%3s  0  0  0  0  0  0  0  0  0  0  0  0" % (
            structure[atom][0],
            structure[atom][1],
            structure[atom][2],
            atom_symbols[atom]))

    for edge in G.edges:
        lines.append("%3s%3s%3s  0" % (new_indexing.index(edge[0]) + 1,
                                       new_indexing.index(edge[1]) + 1,
                                       1))
    lines.append("M  END\n")

    with open(fname, "w") as f:
        f.write("\n".join(lines))


def preprocess_kwargs(kwargs):
    default_keys = {'ignore_elements': [], 'mult': 1.0, 'add_bonds': [], 'remove_bonds': [], 'sdf_name': None}
    for key in default_keys:
        if key not in kwargs:
            kwargs[key] = default_keys[key]


def extract_subgraph(mygraph, kwargs):
    if len(kwargs['ignore_elements']) > 0 and all([isinstance(i, str) for i in kwargs['ignore_elements']]):
        ignore_mode = 'Elements'
    elif len(kwargs['ignore_elements']) > 0 and all([isinstance(i, int) for i in kwargs['ignore_elements']]):
        ignore_mode = 'Indices'
    else:
        assert len(kwargs['ignore_elements']) == 0, "Unexpected contents of ignore_elements: " + repr(kwargs['ignore_elements'])
        ignore_mode = 'None'

    if ignore_mode == 'Elements':
        save_atoms = []
        for node in mygraph.nodes:
            if 'HCarbon' in kwargs['ignore_elements'] and mygraph.nodes[node]['symbol'] == 'H':
                nb_list = list(mygraph.neighbors(node))
                if len(nb_list) == 1 and mygraph.nodes[nb_list[0]]['symbol'] == 'C':
                    continue
            if mygraph.nodes[node]['symbol'] not in kwargs['ignore_elements']:
                save_atoms.append(node)
    elif ignore_mode == 'Indices':
        save_atoms = [node for node in mygraph.nodes if node not in kwargs['ignore_elements' ]]
    else: # if ignore_mode == 'None'
        save_atoms = [node for node in mygraph.nodes]

    for bond in kwargs['add_bonds']:
        atA, atB = bond[0] - 1, bond[1] - 1
        assert atA in save_atoms and atB in save_atoms 
        mygraph.add_edge(atA, atB)
    for bond in kwargs['remove_bonds']:
        atA, atB = bond[0] - 1, bond[1] - 1
        assert atA in save_atoms and atB in save_atoms
        mygraph.remove_edge(atA, atB)

    mysubgraph = mygraph.subgraph(save_atoms)
    return mysubgraph

def generate_connectivity(p, index, kwargs):
    preprocess_kwargs(kwargs)
    
    mygraph = nx.Graph()
    structure = p[index].xyz
    atom_symbols = p.atom_symbols
    
    for i, symbol in enumerate(atom_symbols):
        # if symbol not in ignore_elements:
        mygraph.add_node(i)
        mygraph.nodes[i]['symbol'] = symbol
    for nodeA in range(len(atom_symbols)):
        for nodeB in range(nodeA):
            max_dist = kwargs['mult'] * (RADII[atom_symbols[nodeA]] + RADII[atom_symbols[nodeB]])
            if norm(structure[nodeA] - structure[nodeB]) < max_dist:
                mygraph.add_edge(nodeA, nodeB)
    
    mysubgraph = extract_subgraph(mygraph, kwargs)

    if kwargs['sdf_name'] is not None:
        write_sdf(mysubgraph, p, index, kwargs['sdf_name'])
    return mysubgraph


def generate_connectivity_from_graph(p, molgr, kwargs):
    preprocess_kwargs(kwargs)

    mygraph = nx.Graph()
    atom_symbols = p.atom_symbols
    
    for node, symbol in zip(range(molgr.number_of_nodes()), atom_symbols):
        mygraph.add_node(node)
        mygraph.nodes[node]['symbol'] = symbol
    mygraph.add_edges_from(molgr.edges())

    mysubgraph = extract_subgraph(mygraph, kwargs)
    return mysubgraph


def same_element(n1_attrib, n2_attrib):
    return n1_attrib['symbol'] == n2_attrib['symbol']


def generate_isomorphisms(mysubgraph):
    GM = isomorphism.GraphMatcher(mysubgraph, mysubgraph, node_match=same_element)
    isomorphisms = []
    atoms = list(mysubgraph.nodes)
    for isom in GM.isomorphisms_iter():
        isomorphisms.append([isom[i] for i in atoms])
    return isomorphisms
