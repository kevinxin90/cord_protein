import json
import os
from collections import defaultdict

# mapping between prefix and semantic type
EDGE_TYPE_MAPPING = {
    'CHEBI': 'ChemicalSubstance',
    'CL': 'Cell',
    'DOID': 'DiseaseOrPhenotypicFeature',
    'HGNC': 'Gene',
    'MOP': 'MolecularActivity',
    'NCBITaxon': 'OrganismTaxon',
    'PR': "Protein",
    "SO": "GenomicEntity",
    "UBERON": "AnatomicalEntity"
}

def get_node_type(node, go_mapping):
    """Retrieve the node type based on id.
    
    :param: node: the node id, e.g. GO:00012
    :param: go_mapping: the mapping file between go ID and their semantic type
    """
    prefix = node.split(':')[0]
    if prefix in EDGE_TYPE_MAPPING:
        return EDGE_TYPE_MAPPING[prefix]
    elif node in go_mapping:
        return go_mapping[node]
    return None

def load_cord(data, semantic_type, prefix, go_mapping):
    res = {}
    for edge in data['edges']:
        node1 = edge['node1']['id']
        node2 = edge['node2']['id']
        if node1.split(':')[0] == prefix:
            node1_id = node1.split(':')[1] if prefix == 'HGNC' else node1
            if node1_id not in res:
                res[node1_id] = {
                    prefix.lower(): node1_id,
                    "associated_with": [],
                    "@type": semantic_type
                }
            node_type = get_node_type(node2, go_mapping)
            if node_type:
                node2_prefix = node2.split(':')[0]
                res[node1_id]["associated_with"].append({
                    "@type": node_type,
                    "pmc": edge['evidence'],
                    node2_prefix.lower(): node2 if node_type != "Gene" else node2.split(':')[1]
                })
        if node2.split(':')[0] == prefix:
            node2_id = node2.split(':')[1] if prefix == 'HGNC' else node2
            if node2_id not in res:
                res[node2_id] = {
                    prefix.lower(): node2_id,
                    "associated_with": [],
                    "@type": semantic_type
                }
            node_type = get_node_type(node1, go_mapping)
            if node_type:
                node1_prefix = node1.split(':')[0]
                res[node2_id]["associated_with"].append({
                    "@type": node_type,
                    "pmc": edge['evidence'],
                    node1_prefix.lower(): node1 if node_type != "Gene" else node1.split(':')[1]
                })
    return res

def load_data(data_folder):
    kg_file = os.path.join(data_folder, 'kg.json')
    go_mapping = os.path.join(data_folder, 'go_mapping.json')

    with open(kg_file) as f1:
        kg_data = json.load(f1)

    with open(go_mapping) as f2:
        go_mapping_data = json.load(f2)

    res = load_cord(kg_data, 'Protein', 'PR', go_mapping_data)
    for k, v in res.items():
        v.update({'_id': k})
        yield v
    