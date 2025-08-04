#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 14 09:08:00 2025

@author: guohan
"""

import os
import pandas as pd
from ruamel.yaml import YAML
import string
from tqdm import tqdm


def generate_boltz_yaml(output_dir, polymers, ligand_id, ligand_smiles, affinity=True, constraints=None):
    """
    Generate a single yaml file for Boltz-2.
    :param output_dir: str, path of the output directory.
    :param polymers: list, list of polymers, using the format:
    [{"protein": {"id": "A", "sequence": "KLVFFKLV..."}},
    {"dna":     {"id": "B", "sequence": "ATGCATGC..."}},
    {"rna":     {"id": "C", "sequence": "AUGCGAU..."}}]
    :param ligand_id: str, id of the ligand.
    :param ligand_smiles: str, SMILES of the ligand.
    :param affinity: bool, whether to include affinity in the yaml file.
    :param constraints: list, list of constraints, using the format:
    [{"pocket": {"binder": "X", "contacts": "[ [ A1, 829 ], [ A1, 138 ] ]"}}
    ...]
    """
    os.makedirs(output_dir, exist_ok=True)

    # ligand CHAIN_ID
    polymer_ids = {info["id"] for poly in polymers for info in poly.values()}
    # Select a new unused ligand chain ID (single uppercase letter A-Z)
    for candidate in string.ascii_uppercase:
        if candidate not in polymer_ids:
            ligand_chain_id = candidate
            break
    else:
        raise ValueError("No available CHAIN_ID for ligand")

    # append ligand
    sequences = polymers + [{"ligand": {"id": ligand_chain_id, "smiles": ligand_smiles}}]

    # generate yaml
    data = {
        "version": 1,
        "sequences": sequences
    }
    if affinity:
        data["properties"] = [{"affinity": {"binder": ligand_chain_id}}]
    # define pocket
    if constraints is not None: # [{"pocket": {"binder": "X", "contacts": "[ [ A1, 829 ], [ A1, 138 ] ]"}}]
        for const in constraints:
            for const_type, info in const.items():
                if const_type == 'pocket':
                    info['binder'] = ligand_chain_id
                    # info['contacts'] = str(info['contacts'])
        data["constraints"] = constraints

    output_file = os.path.join(output_dir, f"{ligand_id}.yaml")
    with open(output_file, "w") as f:
        yaml = YAML()
        yaml.dump(data, f)


def batch_generate_boltz_yaml(output_dir, input_file_polymers, input_file_SMILES, id_column_name='ID', smiles_column_name='Cleaned_SMILES', affinity=True):
    """
    Generate a directory of yaml files for Boltz-2.
    :param output_dir: str, path of the output directory.
    :param input_file_polymers: str, path of the input yaml file for polymers.
    :param input_file_SMILES: str, path of the input csv file for SMILES.
    :param id_column_name: str, name of the ID column.
    :param smiles_column_name: str, name of the SMILES column.
    :param affinity: bool, whether to include affinity in the yaml file.
    """
    os.makedirs(output_dir, exist_ok=True)

    # process polymers
    yaml = YAML()
    with open(input_file_polymers, "r") as f:
        doc = yaml.load(f)
    polymers = []
    for entry in doc.get("sequences", []):
        # entry is a dict with a single key in ["protein", "dna", or "rna"], {<ENTITY_TYPE>: {"id": <CHAIN_ID>, "sequence": <sequence>}}
        for entry_type, info in entry.items():
            entry_type = entry_type.lower()
            if entry_type in {"protein", "dna", "rna"}:
                chain_id = info.get("id")
                sequence = info.get("sequence").strip()
                if not chain_id or not sequence:
                    raise ValueError(f"Missing chain_id or sequence in {entry_type} entry")
                polymers.append({entry_type: info})

    # process SMILES
    df = pd.read_csv(input_file_SMILES)
    ligand_IDs = df[id_column_name].tolist()
    ligand_SMILES = df[smiles_column_name].tolist()
    assert len(ligand_IDs) == len(ligand_SMILES), 'Error: Ligand ids and SMILES do not match.'

    # batch generation
    for idx in tqdm(range(len(ligand_IDs))):
        ligand_id, ligand_smiles = ligand_IDs[idx], ligand_SMILES[idx]
        generate_boltz_yaml(output_dir, polymers, ligand_id, ligand_smiles, affinity, constraints=doc.get("constraints", None))


if __name__ == "__main__":
    output_dir = 'tests/test'
    input_file_polymers = 'tests/test_prepare_boltz_yaml_protein.yaml'
    input_file_SMILES = 'tests/test_prepare_boltz_yaml_ligands.csv'
    id_column_name = 'ID'
    smiles_column_name = 'SMILES'
    batch_generate_boltz_yaml(output_dir, input_file_polymers, input_file_SMILES, id_column_name, smiles_column_name,
                              affinity=False)
