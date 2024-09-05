"""
Module to handle ProLIF feature metadata
"""

import logging

import MDAnalysis as mda
from rdkit import Chem
from rdkit.Chem import rdmolops

logger = logging.getLogger("FragFeatures")


class PLFeature:
    """
    Process a ProLIF feature's metadata.
    """

    def __init__(
        self,
        universe: mda.core.universe.Universe,
        mol,
        prolif_feature: dict,
        verbose: bool = False,
        verbose_l2: bool = False,
    ):
        """
        unvierse: MDAnalysis universe object
        prolif_feature: ProLIF feature dictionary
        """
        self.universe = universe
        self.mol = mol
        self.prolif_feature = prolif_feature
        self.interaction_type = next(iter(prolif_feature))
        # print(f"Interaction type: {self.interaction_type}")
        self.interaction_metadata = next(iter(prolif_feature.values()))[0]
        # print(f"Interaction metadata: {self.interaction_metadata}\n\n")
        self.verbose = verbose
        self.verbose_l2 = verbose_l2
        self.get_feature_metadata()

    def get_feature_metadata(self):
        """
        Return the feature metadata for a ProLIF feature.
        """
        if self.verbose_l2:
            print("\nGetting ProLIF feature metadata...")

        self.parent_indices = self.interaction_metadata["parent_indices"]
        self.ligand_indices = self.parent_indices["ligand"]
        self.protein_indices = self.parent_indices["protein"]

        self.interaction_distance = self.interaction_metadata["distance"]

        # TODO: Check if there is additional interaction-specific metadata
        try:
            dha_angle = self.interaction_metadata["DHA_angle"]
        except:
            dha_angle = None
        self.dha_angle = dha_angle

        # NOTE: Not clear what these indices refer to - use parent_indices
        self.indices = self.interaction_metadata["indices"]
        self.ligand_indices_ = self.indices["ligand"]
        self.protein_indices_ = self.indices["protein"]

        self.extract_prot_feat_mdata_from_universe()
        self.extract_lig_feat_mdata_from_rdkit()
        self.get_duck_feature_name()

        if self.verbose_l2:
            print("\n-----------------------")
            print("ProLIF feature metadata")
            print("-----------------------\n")
            print(f"Interaction type: {self.interaction_type}")
            print("")
            print(f"Feature name: {self.duck_feature_name}")
            print("")

            if self.multiple_protein_atoms:
                print(f"Multiple protein atoms: {self.multiple_protein_atoms}")
                print(f"Protein atoms: {self.prot_atom_names}")
                print(f"Protin reisudes: {self.prot_atom_resnames}")
                print(f"Protein ChainIDs: {self.prot_atom_chainID}")
            else:
                print(f"Protein atom: {self.prot_atom_name}")
                print(f"Protein residue: {self.prot_atom_resname}")
                print(f"Protein residue number: {self.prot_atom_resid}")
                print(f"Protein ChainID: {self.prot_atom_chainID}")

            print("")
            print(f"Ligand atoms: {self.ligand_atom_symbols}")
            print(f"Ligand atom charges: {self.ligand_atom_charges}")
            print(f"Ligand atom SMARTS: {self.ligand_atom_smarts}")
            print(f"Ligand substructure SMARTS: {self.smarts_substructure}")
            print("")
            print(f"Interaction distance: {self.interaction_distance:.4f} A")

            if self.dha_angle:
                print(f"DHA angle: {self.dha_angle:.4f} deg")

            print("")
            print(f"Protein indices: {self.protein_indices}")
            print(f"Ligand indices: {self.ligand_indices}")
            print("\n-----------------------\n\n")

    def extract_prot_feat_mdata_from_universe(self):
        """
        Return the DUck feature name for a ProLIF feature.
        """
        if self.verbose_l2:
            print("Extracting protein feature metadata...")
        # MDAnalysis selection
        prot_indices_str = " ".join([str(i) for i in self.protein_indices])
        prot_atomgroup = self.universe.select_atoms(f"index {prot_indices_str}")
        self.prot_atomgroup = prot_atomgroup

        # FIXME: This is a hack to get the first atom in the selection
        # FIXME: Maybe include selection logic to obtain specific atoms
        self.prot_idx = self.protein_indices[0]
        self.prot_atom = self.universe.select_atoms(f"index {self.prot_idx}")

        self.prot_atom_chainID = self.prot_atom.chainIDs[0]
        self.prot_atom_resname = self.prot_atom.resnames[0]
        self.prot_atom_resid = int(self.prot_atom.resids[0])
        self.prot_atom_name = self.prot_atom.names[0]

        # In the case of multiple protein atoms
        self.prot_atom_chainIDs = self.prot_atomgroup.chainIDs
        self.prot_atom_resnames = self.prot_atomgroup.resnames
        self.prot_atom_resids = self.prot_atomgroup.resids
        self.prot_atom_names = self.prot_atomgroup.names

        # NOTE: Multiple protein atoms - currently not suited for DUck
        if len(self.protein_indices) > 1:
            self.multiple_protein_atoms = True
            logger.warning("More than one protein atom associated with protein feature")
        else:
            self.multiple_protein_atoms = False

    def extract_lig_feat_mdata_from_rdkit(self):
        """
        Obtain metadata for ligand feature atoms  using RDKit.
        """
        if self.verbose_l2:
            print("Extracting ligand feature metadata...")
        atoms = {}
        for atom_idx in self.ligand_indices:
            atom = {}
            rdkit_atom_obj = self.mol.GetAtomWithIdx(atom_idx)
            atom_symbol = rdkit_atom_obj.GetSymbol()
            atom_charge = rdkit_atom_obj.GetFormalCharge()
            atom_num_explicit_hs = rdkit_atom_obj.GetNumExplicitHs()
            atom_smarts = rdkit_atom_obj.GetSmarts()

            atom["atom_obj"] = rdkit_atom_obj
            atom["atom_symbol"] = atom_symbol
            atom["atom_charge"] = atom_charge
            atom["atom_num_explicit_hs"] = atom_num_explicit_hs
            atom["atom_smarts"] = atom_smarts

            atoms[f"{atom_symbol}{atom_idx}"] = atom

        self.ligand_atoms = atoms

        self.ligand_atom_symbols = [
            atom["atom_symbol"] for atom in self.ligand_atoms.values()
        ]
        self.ligand_atom_objs = [
            atom["atom_obj"] for atom in self.ligand_atoms.values()
        ]
        self.ligand_atom_charges = [
            atom["atom_charge"] for atom in self.ligand_atoms.values()
        ]
        self.ligand_atom_num_explicit_hs = [
            atom["atom_num_explicit_hs"] for atom in self.ligand_atoms.values()
        ]
        self.ligand_atom_smarts = [
            atom["atom_smarts"] for atom in self.ligand_atoms.values()
        ]

        # SMARTS for the selected substructure
        atom_indices = list(self.ligand_indices)
        try:
            submol = rdmolops.PathToSubmol(self.mol, atom_indices)
            sub_smarts = Chem.MolToSmarts(submol)
            self.smarts_substructure = sub_smarts
        except:
            self.smarts_substructure = None

    def get_duck_feature_name(self):
        """
        Return a feature name prepared for DUck.
        """
        if self.verbose_l2:
            print("\nGetting DUck feature name...")
        feature_name = f"{self.prot_atom_chainID}_{self.prot_atom_resname}_{self.prot_atom_resid}_{self.prot_atom_name}"
        self.duck_feature_name = feature_name

        return self.duck_feature_name
