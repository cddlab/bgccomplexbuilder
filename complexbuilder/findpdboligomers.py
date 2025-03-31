#!/usr/bin/env python3
# %%
"""
PDBの全データからoligomeric stateを取得する
例えばPDB ID: 5E8Hのoligomeric stateはdimeric, 2である。
#
_pdbx_struct_assembly.id                   1
_pdbx_struct_assembly.details              author_and_software_defined_assembly
_pdbx_struct_assembly.method_details       PISA
_pdbx_struct_assembly.oligomeric_details   dimeric
_pdbx_struct_assembly.oligomeric_count     2

APIを叩いても
curl -X GET "https://data.rcsb.org/rest/v1/core/assembly/8DGT/1"
{"pdbx_struct_assembly":{"details":"author_and_software_defined_assembly","id":"1","method_details":"PISA","oligomeric_count":2,"oligomeric_details":"dimeric","rcsb_candidate_assembly":"Y","rcsb_details":"author_and_software_defined_assembly"},"pdbx_struct_assembly_gen":[{"assembly_id":"1","asym_id_list":["A","B","C","D"],"oper_expression":"1","ordinal":1}],"pdbx_struct_oper_list":[{"id":"1","matrix11":1.0,"matrix12":0.0,"matrix13":0.0,"matrix21":0.0,"matrix22":1.0,"matrix23":0.0,"matrix31":0.0,"matrix32":0.0,"matrix33":1.0,"name":"1_555","symmetry_operation":"x,y,z","type":"identity operation","vector1":0.0,"vector2":0.0,"vector3":0.0}],"rcsb_assembly_container_identifiers":{"assembly_id":"1","entry_id":"5E8H","rcsb_id":"5E8H-1","interface_ids":["1"]},"rcsb_assembly_info":{"assembly_id":"1","atom_count":4386,"branched_atom_count":0,"branched_entity_count":0,"branched_entity_instance_count":0,"deuterated_water_count":0,"entry_id":"5E8H","hydrogen_atom_count":0,"modeled_polymer_monomer_count":548,"na_polymer_entity_types":"Other","nonpolymer_atom_count":0,"nonpolymer_entity_count":0,"nonpolymer_entity_instance_count":0,"polymer_atom_count":4127,"polymer_composition":"homomeric protein","polymer_entity_count":1,"polymer_entity_count_dna":0,"polymer_entity_count_rna":0,"polymer_entity_count_nucleic_acid":0,"polymer_entity_count_nucleic_acid_hybrid":0,"polymer_entity_count_protein":1,"polymer_entity_instance_count":2,"polymer_entity_instance_count_dna":0,"polymer_entity_instance_count_rna":0,"polymer_entity_instance_count_nucleic_acid":0,"polymer_entity_instance_count_nucleic_acid_hybrid":0,"polymer_entity_instance_count_protein":2,"polymer_monomer_count":660,"selected_polymer_entity_types":"Protein (only)","solvent_atom_count":259,"solvent_entity_count":1,"solvent_entity_instance_count":2,"unmodeled_polymer_monomer_count":112,"num_interfaces":1,"num_interface_entities":1,"num_homomeric_interface_entities":1,"num_heteromeric_interface_entities":0,"num_isologous_interface_entities":1,"num_heterologous_interface_entities":0,"num_protein_interface_entities":1,"num_na_interface_entities":0,"num_prot_na_interface_entities":0,"total_assembly_buried_surface_area":1791.1230820861406,"total_number_interface_residues":85},"rcsb_id":"5E8H-1","rcsb_struct_symmetry":[{"symbol":"C2","type":"Cyclic","stoichiometry":["A2"],"oligomeric_state":"Homo 2-mer","clusters":[{"members":[{"asym_id":"A","pdbx_struct_oper_list_ids":["1"]},{"asym_id":"B","pdbx_struct_oper_list_ids":["1"]}],"avg_rmsd":0.7324391388321126}],"rotation_axes":[{"start":[0.08568042125161313,76.47996670326157,130.27187021241278],"end":[45.82770172049213,101.85623286493657,127.83938007485918],"order":2}],"kind":"Global Symmetry"}],"rcsb_struct_symmetry_provenance_code":"biojava-7.1.1","rcsb_struct_symmetry_lineage":[{"id":"Global Symmetry.Cyclic.C2","name":"C2","depth":3},{"id":"Global Symmetry.Cyclic.C2.Homo 2-mer","name":"Homo 2-mer","depth":4},{"id":"Global Symmetry","name":"Global Symmetry","depth":1},{"id":"Global Symmetry.Cyclic","name":"Cyclic","depth":2}],"rcsb_latest_revision":{"major_revision":1,"minor_revision":2}}%
[YoshitakaM@cadenza 2:12:59(git)-[main]  ] $ curl -X GET "https://data.rcsb.org/rest/v1/core/assembly/8DGT/1"
{"pdbx_struct_assembly":{"details":"author_defined_assembly","id":"1","oligomeric_count":5,"oligomeric_details":"pentameric","rcsb_candidate_assembly":"Y","rcsb_details":"author_defined_assembly"},"pdbx_struct_assembly_auth_evidence":[{"assembly_id":"1","experimental_support":"gel filtration","id":"1"}],"pdbx_struct_assembly_gen":[{"assembly_id":"1","asym_id_list":["A","B","C","D","E","F","G","H","I","J","K","L","M"],"oper_expression":"1","ordinal":1}],"pdbx_struct_oper_list":[{"id":"1","matrix11":1.0,"matrix12":0.0,"matrix13":0.0,"matrix21":0.0,"matrix22":1.0,"matrix23":0.0,"matrix31":0.0,"matrix32":0.0,"matrix33":1.0,"name":"1_555","type":"identity operation","vector1":0.0,"vector2":0.0,"vector3":0.0}],"rcsb_assembly_container_identifiers":{"assembly_id":"1","entry_id":"8DGT","rcsb_id":"8DGT-1","interface_ids":["1","2","3","4","5","6"]},"rcsb_assembly_info":{"assembly_id":"1","atom_count":11096,"branched_atom_count":0,"branched_entity_count":0,"branched_entity_instance_count":0,"deuterated_water_count":0,"entry_id":"8DGT","hydrogen_atom_count":0,"modeled_polymer_monomer_count":1367,"na_polymer_entity_types":"Other","nonpolymer_atom_count":123,"nonpolymer_entity_count":5,"nonpolymer_entity_instance_count":8,"polymer_atom_count":10973,"polymer_composition":"heteromeric protein","polymer_entity_count":4,"polymer_entity_count_dna":0,"polymer_entity_count_rna":0,"polymer_entity_count_nucleic_acid":0,"polymer_entity_count_nucleic_acid_hybrid":0,"polymer_entity_count_protein":4,"polymer_entity_instance_count":5,"polymer_entity_instance_count_dna":0,"polymer_entity_instance_count_rna":0,"polymer_entity_instance_count_nucleic_acid":0,"polymer_entity_instance_count_nucleic_acid_hybrid":0,"polymer_entity_instance_count_protein":5,"polymer_monomer_count":1904,"selected_polymer_entity_types":"Protein (only)","solvent_atom_count":0,"solvent_entity_count":0,"solvent_entity_instance_count":0,"unmodeled_polymer_monomer_count":537,"num_interfaces":6,"num_interface_entities":6,"num_homomeric_interface_entities":1,"num_heteromeric_interface_entities":5,"num_isologous_interface_entities":1,"num_heterologous_interface_entities":5,"num_protein_interface_entities":6,"num_na_interface_entities":0,"num_prot_na_interface_entities":0,"total_assembly_buried_surface_area":6790.2332833836745,"total_number_interface_residues":384},"rcsb_id":"8DGT-1","rcsb_struct_symmetry":[{"symbol":"C1","type":"Asymmetric","stoichiometry":["A2","B1","C1","D1"],"oligomeric_state":"Hetero 5-mer","clusters":[{"members":[{"asym_id":"B","pdbx_struct_oper_list_ids":["1"]}]},{"members":[{"asym_id":"A","pdbx_struct_oper_list_ids":["1"]}]},{"members":[{"asym_id":"E","pdbx_struct_oper_list_ids":["1"]}]},{"members":[{"asym_id":"C","pdbx_struct_oper_list_ids":["1"]},{"asym_id":"D","pdbx_struct_oper_list_ids":["1"]}],"avg_rmsd":0.7780114698270278}],"kind":"Global Symmetry"},{"symbol":"C1","type":"Asymmetric","stoichiometry":["A2","B2","C1"],"oligomeric_state":"Hetero 5-mer","clusters":[{"members":[{"asym_id":"C","pdbx_struct_oper_list_ids":["1"]},{"asym_id":"D","pdbx_struct_oper_list_ids":["1"]}],"avg_rmsd":0.7664458515080437},{"members":[{"asym_id":"A","pdbx_struct_oper_list_ids":["1"]},{"asym_id":"B","pdbx_struct_oper_list_ids":["1"]}],"avg_rmsd":2.8196605987513372},{"members":[{"asym_id":"E","pdbx_struct_oper_list_ids":["1"]}]}],"kind":"Pseudo Symmetry"},{"symbol":"C2","type":"Cyclic","stoichiometry":["A2"],"oligomeric_state":"Homo 2-mer","clusters":[{"members":[{"asym_id":"C","pdbx_struct_oper_list_ids":["1"]},{"asym_id":"D","pdbx_struct_oper_list_ids":["1"]}],"avg_rmsd":0.7780114698270278}],"rotation_axes":[{"start":[127.45715423048098,136.19772223399838,128.25208623143493],"end":[139.4197062680658,127.72898677653059,82.9092428316647],"order":2}],"kind":"Local Symmetry"}],"rcsb_struct_symmetry_provenance_code":"biojava-7.1.1","rcsb_struct_symmetry_lineage":[{"id":"Pseudo Symmetry.Asymmetric.C1","name":"C1","depth":3},{"id":"Local Symmetry","name":"Local Symmetry","depth":1},{"id":"Local Symmetry.Cyclic.C2.Homo 2-mer","name":"Homo 2-mer","depth":4},{"id":"Global Symmetry.Asymmetric.C1.Hetero 5-mer","name":"Hetero 5-mer","depth":4},{"id":"Pseudo Symmetry","name":"Pseudo Symmetry","depth":1},{"id":"Global Symmetry","name":"Global Symmetry","depth":1},{"id":"Global Symmetry.Asymmetric.C1","name":"C1","depth":3},{"id":"Pseudo Symmetry.Asymmetric","name":"Asymmetric","depth":2},{"id":"Global Symmetry.Asymmetric","name":"Asymmetric","depth":2},{"id":"Pseudo Symmetry.Asymmetric.C1.Hetero 5-mer","name":"Hetero 5-mer","depth":4},{"id":"Local Symmetry.Cyclic","name":"Cyclic","depth":2},{"id":"Local Symmetry.Cyclic.C2","name":"C2","depth":3}],"rcsb_latest_revision":{"major_revision":1,"minor_revision":2}}
"""


import json

import requests


def fetch_pdb_assembly(pdb_id, assembly_id="1"):
    """指定したPDB IDのassembly情報を取得し、整形する"""
    url = f"https://data.rcsb.org/rest/v1/core/assembly/{pdb_id}/{assembly_id}"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        return format_pdb_data(data)
    else:
        print(f"Failed to fetch data for {pdb_id}")
        return None


def format_pdb_data(data):
    """PDB JSONデータを整形して辞書形式で返す"""
    formatted_data = {
        "PDB ID": data.get("rcsb_assembly_container_identifiers", {}).get("entry_id"),
        "Assembly ID": data.get("pdbx_struct_assembly", {}).get("id"),
        "Oligomeric Count": data.get("pdbx_struct_assembly", {}).get(
            "oligomeric_count"
        ),
        "Oligomeric Details": data.get("pdbx_struct_assembly", {}).get(
            "oligomeric_details"
        ),
        "Polymer Composition": data.get("rcsb_assembly_info", {}).get(
            "polymer_composition"
        ),
        "Polymer Entity Count": data.get("rcsb_assembly_info", {}).get(
            "polymer_entity_count"
        ),
        "Polymer Instance Count": data.get("rcsb_assembly_info", {}).get(
            "polymer_entity_instance_count"
        ),
        "Symmetry Type": (
            [sym.get("type") for sym in data.get("rcsb_struct_symmetry", [])]
            if "rcsb_struct_symmetry" in data
            else []
        ),
        "Symmetry Stoichiometry": (
            [sym.get("stoichiometry") for sym in data.get("rcsb_struct_symmetry", [])]
            if "rcsb_struct_symmetry" in data
            else []
        ),
    }

    return formatted_data


# 実行例
pdb_id_list = ["5E8H", "8DGT"]
for pdb_id in pdb_id_list:
    pdb_info = fetch_pdb_assembly(pdb_id)
    if pdb_info:
        print(json.dumps(pdb_info, indent=4))
# %%
