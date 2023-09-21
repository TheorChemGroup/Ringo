import os, sys, platform, json, types
import ringo.pyutils.pyxyz.moltopology
import ringo.pyutils.pyxyz.utils

if platform.system() == "Windows":
    mypath = os.path.dirname(os.path.realpath(__file__))
    if mypath not in sys.path:
        sys.path.insert(0, mypath)
    os.add_dll_directory(mypath)

from .cpppart import cpppart as base
if base.use_pyxyz:
    from .cpppart.cpppart import Confpool
from .pyutils import pyutils
import networkx as nx

DEG2RAD = 0.0174532925199432957692
RAD2DEG = 1 / DEG2RAD

# Factory of Molecule objects
def Molecule(sdf=None, graph_data=None, request_free=None, require_best_sequence=False, clear_feed=True):
    assert (sdf is not None) or (graph_data is not None), "Either 'sdf' or 'graph_data' must be passed as keyword-args"
    assert (sdf is None) or (graph_data is None), "'sdf' or 'graph_data' should not be specified simultaneously"
    
    # Clear feed for new molecule
    if clear_feed:
        base.clear_status_feed()

    m = base.Molecule()

    if request_free is not None:
        request_free = [(i-1, j-1) for i, j in request_free]
        m.set_requested_dofs(request_free)
    
    if require_best_sequence:
        m.require_best_sequence()

    if sdf is not None:
        assert isinstance(sdf, str), "sdf must be an str object (name of your SDF)"
        m.sdf_constructor(sdf, nx, pyutils)
    else: # elif graph_data is not None
        assert isinstance(graph_data, dict), "graph_data must be a dict (see README for details)"
        assert set(graph_data.keys()) == {'graph', 'fixed_dihedrals'}, "graph_data keys must be precisely 'graph' and 'fixed_dihedrals'"
        m.graph_constructor(graph_data, nx, pyutils)
    return m


def get_one_seed():
    return int.from_bytes(os.urandom(3), byteorder="big")

def create_seed_list(size):
    unique_set = set()
    result = []
    while len(result) < size:
        list_element = get_one_seed()
        if list_element not in unique_set:
            unique_set.add(list_element)
            result.append(list_element)
    return result

def cleanup():
    seqcache_path = os.path.join(os.getcwd(), pyutils.CACHE_FILE)
    if os.path.isfile(seqcache_path):
        os.remove(seqcache_path)
        # print(f"Removing cache file {seqcache_path}")

# Work with vdw radii controls
build_flags = base.build_flags
if base.use_overlap_detection:
    get_vdw_radii = base.get_vdw_radii
    set_vdw_radii = base.set_vdw_radii
    set_radius_multiplier = base.set_radius_multiplier

# Work with Ringo status feed
WARNING_CODES = base.warning_codes
for warning_code, warning_line in base.warning_codes.items():
    globals()[warning_code] = warning_line # Declares str variables IK_NOT_APPLIED, SUBOPTIMAL_SOLN_SEQ, UNMET_DOF_REQUEST, etc.
clear_status_feed = base.clear_status_feed
def get_status_feed(important_only=True):
    json_data = base.get_status_feed()
    parsed_data = [
        json.loads(item)
        for item in json_data
    ]
    for item in parsed_data:
        item['important'] = '[important]' in item['subject']
        item['subject'] = item['subject'].replace('[important]', '')
        if len(item['atoms']) > 0:
            item['atoms'] = sorted([idx+1 for idx in item['atoms']])
        else:
            del item['atoms']
    
    if important_only:
        return [item for item in parsed_data if item['important']]
    else:
        return parsed_data

# Summarising statistics for the molecule
def get_molecule_statistics(m):
    graph = m.molgraph_access()
    symbols = m.get_symbols()

    composition = set(
        symbols[atom]
        for atom in graph.nodes
    )
    composition = {element: 0 for element in composition}
    for atom in graph.nodes:
        composition[symbols[atom]] += 1
    
    mcb = [set(c) for c in nx.minimum_cycle_basis(graph)]
    num_cyclic_rotatable_bonds = 0
    num_rotatable_bonds = 0
    for vxA, vxB in graph.edges:
        if len(list(graph.neighbors(vxA))) > 1 and len(list(graph.neighbors(vxB))) > 1:
            num_rotatable_bonds += 1
            for ring_atoms in mcb:
                if vxA in ring_atoms and vxB in ring_atoms:
                    num_cyclic_rotatable_bonds += 1
                    break
    
    dof_list, _ = m.get_ps()

    result = {
        'composition': composition,
        'num_atoms': graph.number_of_nodes(),
        'num_heavy_atoms': len([atom for atom in graph.nodes if symbols[atom] != 'H']),
        'num_bonds': graph.number_of_edges(),
        'num_rotatable_bonds': num_rotatable_bonds,
        'num_cyclic_rotatable_bonds': num_cyclic_rotatable_bonds,
        'largest_macrocycle_size': max([len(item) for item in mcb]),
        'n_flexible_rings': m.get_num_flexible_rings(),
        'n_rigid_rings': m.get_num_rigid_rings(),
        'num_dofs': len(dof_list),
        'cyclomatic_number': graph.number_of_edges() - graph.number_of_nodes() + 1,
    }
    return result
