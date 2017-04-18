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


if __name__ == '__main__':
    ONN = MolGraph().generate('ONN', 'smi')
    ONNO = MolGraph().generate('ONNO', 'smi')
    NN = MolGraph().generate('N#N', 'smi')
    H2O = MolGraph().generate('[OH2]', 'smi')
    NH2OH = MolGraph().generate('NO', 'smi')

    # TEST reaction equivalency
    #all_rxns = scissions(ONNO)
    #new_rxns = scissions(ONNO,all_rxns)
    # assert new_rxns == all_rxns # make sure that reactions don't duplicate

    species = []
    rxns = recursive_scissions(
        [ONNO, ONN, NN, NH3, H2O, NH2OH], finished=species)

    for rxn in rxns:
        print rxn
    print len(rxns)

    for sp in species:
        print sp
    print len(species)
