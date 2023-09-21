import copy, math, json, os, inspect
import networkx as nx
from ..ringo_base import Problem#, log, log_rstart, log_rend
# from ringo_base import Problem#, log, log_rstart, log_rend

CACHE_FILE = 'cached_buildseqs.json'


## GENERATE IG
def calc_ndof(G):
    # addconstr = 0
    ncycles = nx.number_of_edges(G) - nx.number_of_nodes(G) + 1
    nnodes = G.number_of_nodes()
    return nnodes - 3 - 3 * ncycles # - addconstr

def generate_ig(G):
    first_mcb = [set(c) for c in nx.minimum_cycle_basis(G)]
    templist = []
    for item in first_mcb:
        templist.append(sorted(list(item)))
    templist.sort()
    first_mcb = [set(i) for i in templist]

    tempG = nx.Graph()
    for i in range(len(first_mcb)):
        tempG.add_node(i)
        tempG.nodes[i]['atomset'] = first_mcb[i]
        if len(first_mcb[i]) <= 5:
            tempG.nodes[i]['ndof'] = -1 # No DOFs a priory
        else:
            tempG.nodes[i]['ndof'] = len(first_mcb[i]) - 6
    for i in range(tempG.number_of_nodes()):
        for j in range(i + 1, tempG.number_of_nodes()):
            if len(tempG.nodes[i]['atomset'].intersection(tempG.nodes[j]['atomset'])) > 0:
                tempG.add_edge(i, j)
        
    # Unite all a priory rigid rings
    clean_run = False
    while not clean_run:
        clean_run = True
        start_tempG_nodes = list(tempG.nodes)
        for node in start_tempG_nodes:
            if tempG.nodes[node]['ndof'] != -1:
                continue
            # If the part is rigid (ndof == -1)
            for nbnode in list(tempG.neighbors(node)):
                spiro_link = (len(tempG.nodes[node]['atomset']
                                  .intersection(tempG.nodes[nbnode]['atomset'])) == 1)
                lack_dofs = (len(tempG.nodes[nbnode]['atomset']) - 6 - (len(tempG.nodes[node]['atomset'].
                                    intersection(tempG.nodes[nbnode]['atomset'])) - 1) < 0)
                
                if not spiro_link and lack_dofs:
                    # Should unite 'node' with 'nbnode'
                    tempG.nodes[node]['atomset'].update(tempG.nodes[nbnode]['atomset'])
                    for other_nb in tempG.neighbors(nbnode):
                        if other_nb != node:
                            tempG.add_edge(other_nb, node)
                    tempG.remove_node(nbnode)
                    clean_run = False
            if not clean_run:
                break

    mcb = []
    for node in tempG.nodes():
        mcb.append(tempG.nodes[node]['atomset'])
    IG = nx.Graph()
    for i in range(len(mcb)):
        for j in range(i + 1, len(mcb)):
            common_atoms = mcb[i].intersection(mcb[j])
            if len(common_atoms) > 0:
                IG.add_edge(i, j)
                IG[i][j]['link'] = list(common_atoms)
    if IG.number_of_edges() == 0 and len(mcb) == 1:
        IG.add_node(0)
    assert IG.number_of_nodes() != 0, "Cycle graph has 0 nodes"
    
    for i in range(len(mcb)):
        fraggraph = nx.Graph()
        for edge in list(G.edges):
            if edge[0] in mcb[i] and edge[1] in mcb[i]:
                issingle = (G[edge[0]][edge[1]]['type'] == 1)
                fraggraph.add_edge(*edge, single=issingle)
        IG.nodes[i]['graph'] = fraggraph
        IG.nodes[i]['NDOF'] = calc_ndof(fraggraph)

    for nA, nB in IG.edges:
        comm_bonds = []
        for edge in IG.nodes[nA]['graph'].edges:
            if IG.nodes[nB]['graph'].has_edge(*edge):
                comm_bonds.append(sorted(list(edge)))
        IG[nA][nB]['common_bonds'] = sorted(comm_bonds)
    return IG


## GENERATE BUILD SEQUENCE
def log_problem(mygraph, consbonds, seq):
    log_bondtypes = []
    for edge in mygraph.edges:
        if mygraph[edge[0]][edge[1]]['type']:
            value = 1
        else:
            value = 0
        log_bondtypes.append([*edge, value])
    
    cb_repr=None
    if len(consbonds) == 0:
        cb_repr = repr([])
    elif isinstance(consbonds[0], list):
        cb_repr = repr(sorted(consbonds))
    # if cb_repr is not None:
    #     log_rstart("problem_init")
    #     log(" IN:: nodes = " + repr(sorted(list(mygraph.nodes))))
    #     log(" IN:: consbonds = " + cb_repr)
    #     log(" IN:: seq = " + repr(seq))
    #     log(" CHECK:: edges_SET2 = " + repr(list(mygraph.edges)))
    #     log(" CHECK:: bondtypes_SET2 = " + repr(log_bondtypes))
    #     log_rend("problem_init")

def sync_types(IG, node):
    for nb in IG.neighbors(node):
        link = IG[node][nb]['link']
        if len(link) == 2:
            if not IG.nodes[node]['graph'][link[0]][link[1]]['type']:
                IG.nodes[nb]['graph'][link[0]][link[1]]['type'] = False
        elif len(link) > 2:
            for edge in IG.nodes[node]['graph'].edges():
                if edge[0] in link and edge[1] in link and \
                        not IG.nodes[node]['graph'][edge[0]][edge[1]]['type']:
                    IG.nodes[nb]['graph'][edge[0]][edge[1]]['type'] = False
    for nb in IG.neighbors(node):
        link = IG[node][nb]['link']
        if len(link) == 2:
            if not IG.nodes[node]['graph'][link[0]][link[1]]['unrestrained']:
                IG.nodes[nb]['graph'][link[0]][link[1]]['unrestrained'] = False
        elif len(link) > 2:
            for edge in IG.nodes[node]['graph'].edges():
                if edge[0] in link and edge[1] in link and \
                        not IG.nodes[node]['graph'][edge[0]][edge[1]]['unrestrained']:
                    IG.nodes[nb]['graph'][edge[0]][edge[1]]['unrestrained'] = False

def process_seq(seq, IG, requested_free, final=True):
    found_bonds = [False] * len(requested_free)
    for node in seq:
        G = IG.nodes[node]['graph']
        for edge in G.edges:
            G[edge[0]][edge[1]]['type'] = G[edge[0]][edge[1]]['single']
            G[edge[0]][edge[1]]['unrestrained'] = True
        for i, bond in enumerate(requested_free):
            if G.has_edge(*bond):
                found_bonds[i] = True
    assert False not in found_bonds

    LG = nx.DiGraph(directed=True)
    for node in seq:
        LG.add_node(node)
        for nb in IG.neighbors(node):
            if LG.has_node(nb):
                LG.add_edge(nb, node)

    for node in LG.nodes:
        nb_list = list(LG.neighbors(node))
        LG.nodes[node]['dep_max'] = len(nb_list)
        LG.nodes[node]['dep_cur'] = 0
        LG.nodes[node]['ikprob'] = None

    LG = LG.reverse(copy=False)
    for node in LG.nodes:
        nb_list = list(LG.neighbors(node))
        lbs = []
        for nb in nb_list:
            for bond in IG[node][nb]['common_bonds']:
                if bond not in lbs:# and [bond[1], bond[0]] not in lbs:
                    lbs.append(bond)
        LG.nodes[node]['linking_bonds'] = lbs

    done = False
    onemorerun = True

    reached_ends = 0
    num_of_ends = len([node for node in LG.nodes if len(list(LG.neighbors(node))) == 0])
    while not done:
        for node in seq:
            if LG.nodes[node]['dep_cur'] == LG.nodes[node]['dep_max']:
                if LG.nodes[node]['ikprob'] is None:
                    current_request_free = []
                    for bond in requested_free:
                        if IG.nodes[node]['graph'].has_edge(*bond) and \
                           list(bond) not in LG.nodes[node]['linking_bonds'] and \
                           [bond[1], bond[0]] not in LG.nodes[node]['linking_bonds']:
                            current_request_free.append(bond)
                    LG.nodes[node]['ikprob'] = Problem(IG.nodes[node]['graph'], LG.nodes[node]['linking_bonds'], current_request_free, final)
                    sync_types(IG, node)
                elif not LG.nodes[node]['ikprob'].recheck_method():
                    sync_types(IG, node)
                    onemorerun = True

                if len(list(LG.neighbors(node))) == 0:
                    reached_ends += 1
                
                if reached_ends == num_of_ends:
                    reached_ends = 0
                    done = True
                    if onemorerun:
                        onemorerun = False
                        done = False
                        for node in LG.nodes():
                            LG.nodes[node]['dep_cur'] = 0
                    continue

                for nb in LG.neighbors(node):
                    LG.nodes[nb]['dep_cur'] += 1
                LG.nodes[node]['dep_cur'] += 1 # To avoid repetitions for the same node
    
    if final:
        # We have TLC applicability warnings prepared during iterations of 'while'
        # but they are recored afterwards.
        for node in seq:
            LG.nodes[node]['ikprob'].record_warning()
    
        blocked_graph = nx.DiGraph(directed=True)
        blocked_graph.add_edges_from(LG.edges)
        nonblocked_nodes = []
        for node in blocked_graph.nodes:
            frag_graph = IG.nodes[node]['graph']

            # Not interested in rings that are processed by TLC
            unrestrained = False
            for bondA, bondB in frag_graph.edges:
                if frag_graph[bondA][bondB]['unrestrained']:
                    unrestrained = True
                    break
            assert ((LG.nodes[node]['ikprob'].method == 1) and unrestrained) or \
                ((LG.nodes[node]['ikprob'].method != 1) and not unrestrained)
            if unrestrained:
                nonblocked_nodes.append(node)
                continue
            
            # Not interested in rings that are rigid purely by their NDOFs
            nedges = frag_graph.number_of_edges()
            nnodes = frag_graph.number_of_nodes()
            nrings = nedges - nnodes + 1
            nfixed = 0
            for edge in frag_graph.edges:
                if not frag_graph[edge[0]][edge[1]]['single']:
                    nfixed += 1
            if nedges - nrings * 6 - nfixed < 0:
                # print(f"Node #{node} is rigid purely by NDOFs. nedges={nedges} nnodes={nnodes} nrings={nrings} nfixed={nfixed}")
                nonblocked_nodes.append(node)
        
        # Spiro-linkages do not create any dependencies between rings
        # Removing them will make generated warnings clearer
        remove_edges = []
        for nodeA, nodeB in blocked_graph.edges:
            if len(IG[nodeA][nodeB]['common_bonds']) == 0:
                remove_edges.append((nodeA, nodeB))
        blocked_graph.remove_edges_from(remove_edges)

        # Remove the fragments that are totally not flexible or IK will be applied
        blocked_graph.remove_nodes_from(nonblocked_nodes)

        if blocked_graph.number_of_nodes() > 0:
            blocked_undirected = nx.Graph()
            blocked_undirected.add_edges_from(blocked_graph.edges)

            # Each connected component represents a potentially flexible system
            # that falls under limitations for our polycycle processing method
            for component_nodes in nx.connected_components(blocked_undirected):
                component_graph = IG.subgraph(component_nodes)
                component_moledges = []
                for node in component_graph.nodes:
                    component_moledges += IG.nodes[node]['graph'].edges

                component_molgraph = nx.Graph()
                component_molgraph.add_edges_from(component_moledges)
                for node in component_graph.nodes:
                    for edge in IG.nodes[node]['graph'].edges:
                        component_molgraph[edge[0]][edge[1]]['single'] = IG.nodes[node]['graph'][edge[0]][edge[1]]['single']

                nedges = component_molgraph.number_of_edges()
                nnodes = component_molgraph.number_of_nodes()
                nrings = nedges - nnodes + 1
                nfixed = 0
                for edge in component_molgraph.edges:
                    if not component_molgraph[edge[0]][edge[1]]['single']:
                        nfixed += 1
                
                if nedges - nrings * 6 - nfixed >= 0:
                    message_data = {
                        "message": f"Unable to sample kinematically flexible {len(component_nodes)}-cyclic fragment ({component_molgraph.number_of_nodes()} atoms in the fragment). This fragment will be kept as it is in the starting conformation.",
                        "subject": Problem.get_warning_codes()["IK_NOT_APPLIED"]+"[important]",
                        "atoms": list(component_molgraph.nodes),
                        "file": __file__,
                        "line": inspect.currentframe().f_lineno,
                    }
                    Problem.add_message_to_feed(json.dumps(message_data))
                else:
                    print(f"[CRITICAL ERROR] This is certainly a bug. Encountered kinematically rigid polycycle where should have not. {__file__}:{inspect.currentframe().f_lineno}")
        
        # Add warnings if there are unfulfilled DOF requests
        for node in LG.nodes:
            unfulfilled = LG.nodes[node]['ikprob'].get_unfulfilled_requests()
            for bond in unfulfilled:
                message_data = {
                        "message": f"Unable to enforce the bond {[i+1 for i in bond]} to be a free DOF of the ring(s) {repr([i+1 for i in IG.nodes[node]['graph'].nodes])}",
                        "subject": Problem.get_warning_codes()["UNMET_DOF_REQUEST"]+"[important]",
                        "atoms": list(bond),
                        "file": __file__,
                        "line": inspect.currentframe().f_lineno,
                    }
                Problem.add_message_to_feed(json.dumps(message_data))
    return LG

def strats_with_headnode(headnode, IG, requested_free):
    done_seqs = []
    prospects = [[headnode]]
    while len(prospects) > 0:
        cur_prospect = prospects.pop()
        # TODO Check cur_prospect is okay

        if len(cur_prospect) < IG.number_of_nodes():
            for node in IG.nodes:
                if node not in cur_prospect:
                    prospects.append(cur_prospect + [node])
        elif cur_prospect not in done_seqs:
            done_seqs.append(cur_prospect)

    scores = []
    # for node in IG.nodes:
    #     print(f"Nodes of CP#{node} = {repr(list(IG.nodes[node]['graph'].nodes))}")
    # print(f"SEQS = {repr(done_seqs)}")
    for seq_idx, seq in enumerate(done_seqs):
        # print(f"----------- SEQ #{seq_idx} -----------")
        LG = process_seq(seq, IG, requested_free, final=False)

        score = 0.0
        for node in LG.nodes:
            if LG.nodes[node]['ikprob'].method != 0:
                score += 1.0
            unfulfilled = LG.nodes[node]['ikprob'].get_unfulfilled_requests()
            score -= len(unfulfilled)*0.01
        scores.append(score)
    
    best_score = max(scores)
    best_strats = [done_seqs[i] for i in range(len(scores)) if scores[i] == best_score]

    log_atoms = set()
    for node in IG.nodes:
        G = IG.nodes[node]['graph']
        log_atoms.update(set(G.nodes))
    # log_rstart("stratheadnode_final")
    # log(" IN:: nodes = " + repr(sorted(list(IG.nodes))))
    # log(" IN:: headnode = " + repr(headnode))
    # log(" IN:: atoms = " + repr(sorted(list(log_atoms))))
    # log(" CHECK:: headatoms = " + repr(sorted(list(IG.nodes[headnode]['graph'].nodes))))
    # log(" CHECK:: bestscore = " + repr([best_score]))
    # log(" CHECK:: best_strats = " + repr(best_strats))
    # log_rend("stratheadnode_final")

    return best_strats, best_score

def get_strategy_simple(IG):
    LG = nx.Graph()
    max_dof_value = None
    max_dof_node = None
    for node in IG.nodes:
        LG.add_node(node)
        if max_dof_value is None or max_dof_value < IG.nodes[node]['NDOF']:
            max_dof_value = IG.nodes[node]['NDOF']
            max_dof_node = node
    
    done = False
    onemorerun = True
    seq = [max_dof_node] # Built in reverse order
    while not done:
        done = True
        for node in seq:
            for nb in IG.neighbors(node):
                if nb not in seq:
                    LG.add_edge(node, nb)
                    seq.append(nb)
                    onemorerun = True
                    done = False
        if onemorerun:
            done = False
            onemorerun = False
    seq = [i for i in reversed(seq)] # Assemble in normal order
    return seq

def get_strategy(IG, filename, requested_free, require_best_sequence):
    # Don't waste time in the most trivial cases
    if IG.number_of_nodes() == 1:
        assert IG.has_node(0)
        return [0]

    cached_data = {}
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE, 'r') as f:
            cached_data = json.load(f)
    
    ig_data = {
        'ignodes': sorted(list(IG.nodes)),
        'igedges': sorted([sorted(list(edge))for edge in IG.edges]),
        'cycles': [sorted(list(IG.nodes[node]['graph'].nodes)) for node in IG.nodes]
    }
    
    # Check if the current cyclic part is cached in the file
    present = False
    cached_index = None
    if filename in cached_data:
        for cur_index, cur_data in enumerate(cached_data[filename]):
            if cur_data['id'] == ig_data:
                present = True
                cached_index = cur_index
    
    # Get build seq from cached_data
    if present:
        build_seq = cached_data[filename][cached_index]['seq']
        # print(f'Build sequence has been retrieved from {CACHE_FILE}')
        return build_seq

    # If cache is not available, compute build seq from scratch
    # print(f'Number of required iterations = {math.factorial(IG.number_of_nodes())}')
    if not require_best_sequence and math.factorial(IG.number_of_nodes()) > 40000:
        # Post a warning that bruteforce is too demanding
        atoms_set = set()
        for node in IG.nodes:
            G = IG.nodes[node]['graph']
            atoms_set.update(set(G.nodes))
        add_warning = ""
        if len(requested_free) > 0:
            add_warning = " Unable to enforce the requested dihedral(s)."
        message_data = {
            "message": f"Cannot use bruteforce to find the best solution sequence. If there are no further warnings IK will still be applicable.{add_warning} In case if this causes any problems consider using require_best_sequence=True.",
            "subject": Problem.get_warning_codes()["SUBOPTIMAL_SOLN_SEQ"],
            "atoms": sorted(list(atoms_set)),
            "file": __file__,
            "line": inspect.currentframe().f_lineno,
        }
        Problem.add_message_to_feed(json.dumps(message_data))

        res = get_strategy_simple(IG)
    else:
        # If bruteforce is feasible
        best_strat = None
        best_score = None

        for headnode in IG.nodes:
            strats, score = strats_with_headnode(headnode, IG, requested_free)
            if best_score is None or best_score < score:
                best_score = score
                best_strat = strats[0]
            
        log_atoms = set()
        for node in IG.nodes:
            G = IG.nodes[node]['graph']
            log_atoms.update(set(G.nodes))
        res = best_strat
    
    # Cache the build seq
    if filename not in cached_data:
        cached_data[filename] = []    
    cached_data[filename].append({
        'id': ig_data,
        'seq': res,
    })
    with open(CACHE_FILE, 'w') as f:
        json.dump(cached_data, f)
    return res
