from imolecule.notebook import generate
from imolecule.format_converter import convert
import networkx as nx
import networkx.algorithms.isomorphism as iso
from data_structures import MolGraph, RxnGraph


def unique_graphs(graph_list):
    uniques = []
    for ga in graph_list:
        if True not in [ga == gb for gb in uniques]:
            uniques.append(ga)
    return uniques


def scissions(reactant, rxns=[]):

    new_mols = []
    for edge in reactant.edges():
        new_reactant = reactant.copy()
        new_reactant.remove_edge(*edge)
        new_mols.append(new_reactant)

    # check isomorphism to identify other duplicates
    unique_mols = unique_graphs(new_mols)

    # use k-component structure to identify fragments
    all_frags = []
    for mol_a in unique_mols:
        if len(mol_a) == 2:  # bimolecular species can't be decomposed with k_components
            node_a, node_b = mol_a.nodes()
            subg_a = mol_a.subgraph([node_a])
            subg_b = mol_a.subgraph([node_b])
            rxn_list = [[reactant], [subg_a, subg_b]]
        else:
            k = nx.k_components(mol_a)
            all_nodes = set(mol_a.nodes())
            for nk in k:
                rxn = []
                frags = k[nk]
                for f in frags:
                    subg = mol_a.subgraph(f)
                    all_nodes -= f
                    rxn.append(subg)
                    all_frags.append(subg)

                if all_nodes:
                    # create fragment from remaining nodes
                    subg = mol_a.subgraph(all_nodes)
                    rxn.append(subg)

                rxn_list = [[reactant], rxn]

        new_rxn = RxnGraph()
        new_rxn.from_rxn_list([rxn_list])

        rxns = unique_graphs(rxns + [new_rxn])

    return rxns


def recursive_scissions(reactant_list, network=[], finished=[]):

    species = []

    for reactant in reactant_list:
        finished.append(reactant)
        rxns = scissions(reactant)
        for rxn in rxns:
            species += [
                rxn.node[n]['graph'] for n in rxn.nodes() if (
                    rxn.node[n]['type'] == 'molecule' and rxn.node[n]['molecule_type'] == 'product' and len(
                        rxn.node[n]['graph']) > 1)]

        network += rxns
        network = unique_graphs(network)

    species = unique_graphs(species)
    species = [sp for sp in species if sp not in finished]

    if not species:
        return network

    print 'Recursing on:'
    print [str(sp) for sp in species]
    return unique_graphs(
        network +
        recursive_scissions(
            species,
            network=network,
            finished=finished))


def get_mechanisms(rxn_net, starting_species, terminal_species, max_depth):

    def branches(species):
        b = []
        for rxn in rxn_net.successors(species):
            prods = rxn_net.successors(rxn)
            b.append(prods)
        return b

    for species in starting_species:
        print(species)
        for i, b in enumerate(branches(species)):
            print('Branch {}'.format(i))
            print([str(bi) for bi in b])

if __name__ == '__main__':

    ONNO = MolGraph().generate('ONNO','smi') #generate ONNO
    ONN = MolGraph().generate('ONN','smi') #generate ONN
    N2H4 = MolGraph().generate('NN','smi') #generate N2H4
    NH3 = MolGraph().generate('[NH3]','smi') #generate NH3
    H2O = MolGraph().generate('[OH2]','smi') #generate H2O
    O2 = MolGraph().generate('O=O','smi') #generate O2
    N2 = MolGraph().generate('N#N','smi') #generate N2
    H2 = MolGraph().generate('[HH]','smi') #generate H2
    NH2OH = MolGraph().generate('NO','smi') #generate NH2OH
    NO = MolGraph().generate('[N]=O','smi') #generate NO

    global_reactants = [H2O, N2, O2]
    key_intermediates = [ONNO, ONN, N2H4]
    global_products = [NO,NH3, NH2OH]

    rxns = recursive_scissions(
        [ONN,O2])

    rxns_list = []
    for rxn in rxns:
	products = [rxn.node[p]['graph'] for p in rxn.nodes() if rxn.node[p]['type']=='molecule' and rxn.node[p]['molecule_type']=='product']
	reactants = [rxn.node[p]['graph'] for p in rxn.nodes() if rxn.node[p]['type']=='molecule' and rxn.node[p]['molecule_type']=='reactant']
	rxns_list.append([reactants, products])

    all_rxns = RxnGraph()
    all_rxns.from_rxn_list(rxns_list)

    all_rxns.classify_rxns()

    rxns = [n for n in all_rxns.nodes() if all_rxns.node[n]['type'] == 'reaction']

    for rxn in rxns:
	reacts, prods = all_rxns.get_reactants_products(rxn)
	rxn_type = all_rxns.node[rxn]['reaction_type']
	react_composition = {}
	prod_composition = {}
	for r in reacts:
	    react_composition = r.composition(react_composition) #cumulative reactant composition
	for r in prods:
	    prod_composition = r.composition(prod_composition) #cumulative product composition
	    
	if react_composition.get('N',0) == 2:
	    #reactions before N-N bond scission
	    if rxn_type not in ['N-N_scission','H-H_scission']:
		all_rxns.reverse_rxn(rxn,remove=True) #make N2-H and N2-O scissions into couplings
	    elif True in [r in global_products for r in reactants]:
		all_rxns.reverse_rxn(rxn,remove=True) #Make sure global products are also local products

    get_mechanisms(all_rxns, ['NN','OO'], ['NO'], 5)
