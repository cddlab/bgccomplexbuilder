# https://data.rcsb.org/graphql/index.html

{
  entries(entry_ids: ["3WWN"]) {
    rcsb_id
    entry {
      id
    }
    pubmed {
      rcsb_pubmed_container_identifiers {
        pubmed_id
      }
      rcsb_pubmed_central_id
      rcsb_pubmed_doi
      rcsb_pubmed_abstract_text
      rcsb_pubmed_affiliation_info
    }
    assemblies {
      rcsb_assembly_container_identifiers {
        assembly_id
        entry_id
        rcsb_id
        interface_ids
      }
      pdbx_struct_assembly {
        details
        id
        oligomeric_count
        oligomeric_details
        rcsb_candidate_assembly
        rcsb_details
      }
      pdbx_struct_assembly_auth_evidence {
        experimental_support
      }
      rcsb_struct_symmetry {
        kind
        type
        symbol
        oligomeric_state
        stoichiometry
      }
      rcsb_assembly_info {
        assembly_id
        atom_count
        branched_atom_count
        branched_entity_count
        branched_entity_instance_count
        deuterated_water_count
        entry_id
        hydrogen_atom_count
        modeled_polymer_monomer_count
        na_polymer_entity_types
        nonpolymer_atom_count
        nonpolymer_entity_count
        nonpolymer_entity_instance_count
        polymer_atom_count
        polymer_composition
        polymer_entity_count
        polymer_entity_count_DNA
        polymer_entity_count_RNA
        polymer_entity_count_nucleic_acid
        polymer_entity_count_nucleic_acid_hybrid
        polymer_entity_count_protein
        polymer_entity_instance_count
        polymer_entity_instance_count_DNA
        polymer_entity_instance_count_RNA
        polymer_entity_instance_count_nucleic_acid
        polymer_entity_instance_count_nucleic_acid_hybrid
        polymer_entity_instance_count_protein
        polymer_monomer_count
        selected_polymer_entity_types
        solvent_atom_count
        solvent_entity_count
        solvent_entity_instance_count
        unmodeled_polymer_monomer_count
        num_interfaces
        num_interface_entities
        num_homomeric_interface_entities
        num_heteromeric_interface_entities
        num_isologous_interface_entities
        num_heterologous_interface_entities
        num_protein_interface_entities
        num_na_interface_entities
        num_prot_na_interface_entities
        total_assembly_buried_surface_area
        total_number_interface_residues
      }
    }
  }
}