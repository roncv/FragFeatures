"""
Module for generating compound/fragment fingerprints from protein-ligand complexes.
"""

import logging

# TODO: Rename this to ...
import os

import MDAnalysis as mda
import molparse as mp
import numpy as np
import prolif as plf

# from molparse.rdkit.features import FEATURE_FAMILIES, COMPLEMENTARY_FEATURES
from molparse.rdkit import draw_mols

# RDKIT imports
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, Draw, MolFromSmiles

from FragFeatures.core.tools import sanitise_mol
from FragFeatures.fingerprint.prolif_feature import PLFeature
from FragFeatures.sanitisation.protein_preparation import ProteinPreparation
from FragFeatures.utils import timefunction

logging.getLogger("MDAnalysis").setLevel(logging.WARNING)
logger = logging.getLogger("FragFeatures")  # NOTE: Implement this a bit more



# TODO: Evaluate cutoffs and feature selection
# TODO: Consider a an alternative way to define directory structures
class Pose:
    """
    Represents a pose of a compound in a target.
    """

    def __init__(
        self,
        target_dir,
        compound_code,
        selected_interactions="hbonds+",
        verbose=False,
        verbose_l2=False,
    ):
        self.target_dir = target_dir
        self.compound_code = compound_code
        self.id = compound_code  # TODO: for now: what is the id for a pose? Ask Max
        self.protein_path = (
            f"{target_dir}/aligned_files/{compound_code}/{compound_code}_apo.pdb"
        )
        self.mol_path = (
            f"{target_dir}/aligned_files/{compound_code}/{compound_code}_ligand.mol"
        )
        self.smi_path = (
            f"{target_dir}/aligned_files/{compound_code}/{compound_code}_ligand.smi"
        )
        self.verbose = verbose
        self.verbose_l2 = verbose_l2
        self.plf_interaction_types = self.define_prolif_interation_types(selected_interactions)
        self.protein_system = mp.parse(
            self.protein_path, verbosity=False
        ).protein_system

        self.sparse_fingerprint = (
            None  # TODO: Is this the way to store the fingerprint? Ask Max
        )
        self.sparse_fingerprint_ext = None
        self.fingerprint = None
        self.fingerprint_ext = None

        self._mol = None  # TODO: Why is there an underscore here? Ask Max

    def define_prolif_interation_types(self, selected_interactions):
        """
        Define interaction types for the pose from groupings of prolif feature families.
        """
        feature_families = [
            'Anionic',
            'CationPi',
            'Cationic',
            'EdgeToFace',
            'FaceToFace',
            'HBAcceptor',
            'HBDonor',
            'Hydrophobic',
            'MetalAcceptor',
            'MetalDonor',
            'PiCation',
            'PiStacking',
            'VdWContact',
            'XBAcceptor',
            'XBDonor'
            ]

        valid_interactions_selections = [
            "hbonds",
            "hbonds+",
            "-aromatic",
            "all",
            ]

        if selected_interactions not in valid_interactions_selections:
            raise ValueError(
                f"Invalid interaction types: {selected_interactions}. Valid options are: {valid_interactions_selections}"
            )

        if selected_interactions == "hbonds":
            interaction_types = [
                "HBDonor",
                "HBAcceptor",
                "XBAcceptor",
                "XBDonor",
                ]
        elif selected_interactions == "hbonds+":
            interaction_types = [
                "HBDonor",
                "HBAcceptor",
                "XBAcceptor",
                "XBDonor",
                "Anionic",
                "Cationic",
                "PiCation",
                ]
        elif selected_interactions == "-aromatic":
            interaction_types = [
                "HBDonor",
                "HBAcceptor",
                "XBAcceptor",
                "XBDonor",
                "Anionic",
                "Cationic",
                "PiCation",
                "Hydrophobic",
                "MetalAcceptor",
                "MetalDonor",
                "VdWContact",
                ]
        elif selected_interactions == "all":
            interaction_types = feature_families
        
        if self.verbose:
            print(f"\nProLIF Interaction Types: {interaction_types}")

        return plf.Fingerprint(interaction_types)

    @timefunction
    def calculate_prolif_fp(self, output_dir, rdkit_protein=True):
        """
        Calculate the ProLIF fingerprint for a pose.
        """
        if self.verbose:
            print(f"\nCalculating ProLIF fingerprint for {self.compound_code}...")
            print(f"\nProtein path: {self.protein_path}")
            print(f"\nLigand path: {self.mol_path}")

        protein_prep = ProteinPreparation(
            protein_path=self.protein_path,
            output_dir=output_dir,
            protein_id=self.compound_code,
            minimize=True,
            pH=7.8, # NOTE: This is hardcoded
            verbose=self.verbose,
            verbose_l2=self.verbose_l2,
        )
        protein_prep.prepare_protein()
        prepared_protein_path = protein_prep.get_prepared_protein_path()
        protein_filename = protein_prep.get_protein_filename()
        self.protein_termini = protein_prep.termini
        self.prepared_protein_path = prepared_protein_path
        self.protein_filename = protein_filename

        if rdkit_protein:
            try:
                rdkit_prot = Chem.MolFromPDBFile(prepared_protein_path, removeHs=False)
                protein_mol = plf.Molecule(rdkit_prot)
            except:
                logger.warning(
                    f"Could not load protein from {prepared_protein_path} using RDKit. Falling back to MDAnalysis."
                )
                u = mda.Universe(prepared_protein_path)
                protein_mol = plf.Molecule.from_mda(u)
        else:
            u = mda.Universe(prepared_protein_path)
            protein_mol = plf.Molecule.from_mda(u)

        fp = self.plf_interaction_types

        self.plf_fp_obj = fp

        # FIXME: Make this a property
        self.ligand_mol = self.ligand_preparation(output_dir=output_dir)
        ligand = plf.Molecule.from_rdkit(self.ligand_mol)

        fp.run_from_iterable([ligand], protein_mol)  # get fingerprint
        df = fp.to_dataframe()

        if self.verbose:
            print(f"\nProLIF Interactions Dataframe\n\n{df.T}\n")

        # ProLIF only identifies using ligand and protein residues
        res_lig_interactions = []
        for lig_res, prot_res, interaction in df.columns:
            prot_res_result = fp.ifp[0][(lig_res, prot_res)]
            res_lig_interactions.append((lig_res, prot_res))
        # Remove duplicates
        unique_res_lig_interactions = list(set(res_lig_interactions))

        self.number_of_interactions_plf = len(df.columns)\
        # FIXME: How to deal with this better
        if self.number_of_interactions_plf == 0:
            logger.warning(
                f"No interactions found for {self.compound_code} using ProLIF."
            )
        self.plf_fp_keys = unique_res_lig_interactions

        if self.verbose:
            print(f"\nNumber of interactions: {self.number_of_interactions_plf}")
            print(f"\nPLF FP Keys\n\n{self.plf_fp_keys}\n\n")

        interactions = []
        for lig_res, prot_res in self.plf_fp_keys:
            prot_res_result = fp.ifp[0][(lig_res, prot_res)]  # access the interaction
            if len(prot_res_result) > 1:
                for key, value in prot_res_result.items():
                    prot_res_result_i = {}
                    prot_res_result_i[key] = value
                    interactions.append(prot_res_result_i)
            else:
                interactions.append(prot_res_result)

        # This is a list of dictionaries, each being a single interaction
        self.plf_fp = interactions  # TODO: Turn this into a property

        if self.verbose:
            print(
                f"Extracting metadata from {self.number_of_interactions_plf} features..."
            )

        u = mda.Universe(prepared_protein_path)
        duck_feature_names = []
        plf_feature_objs = []
        for feature in self.plf_fp:
            plf_feature = PLFeature(
                universe=u,
                mol=self.ligand_mol,
                prolif_feature=feature,
                verbose=self.verbose,
                verbose_l2=self.verbose_l2,
            )
            duck_feature_name = plf_feature.duck_feature_name

            duck_feature_names.append(duck_feature_name)
            plf_feature_objs.append(plf_feature)

        self.duck_feature_names_plf = duck_feature_names
        self.plf_feature_objs = plf_feature_objs

        if self.verbose:
            print(f"\nDUck feature names:\n{duck_feature_names}")

    def get_interactive_ligand_network(self, fp, output_dir):
        """
        Return the interactive protein-ligand interaction network plot for the pose
        """
        try:
            fp.plot_lignetwork(self.ligand_mol, kind="frame", frame=0)
            network_plot = fp.plot_lignetwork(self.ligand_mol, kind="frame", frame=0)
        except:
            # FIXME: Raise an exception if there are no interactions??
            # logger.warning(
            #     f"Could not plot the interactive protein-ligand interaction network for {self.compound_code}."
            # )
            return

        # Assuming 'output_dir' is your directory path where you want to save the HTML file
        html_file_path = os.path.join(
            output_dir, f"{self.compound_code}_interactive_plot.html"
        )

        with open(html_file_path, "w") as file:
            file.write(network_plot.data)  # Assuming 'data' contains the HTML content

    def ligand_preparation(self, output_dir=None):
        """
        Prepare the ligand for featrure analysis.
        """
        m = Chem.MolFromMolFile(self.mol_path)
        self.smiles = Chem.MolToSmiles(m)

        # m = AllChem.AssignBondOrdersFromTemplate(m, )

        if self.verbose_l2:
            print(f"\n{self.compound_code} NumAtoms: {m.GetNumAtoms()}")
        # Chem.AllChem.EmbedMolecule(m)

        m = Chem.AddHs(m, addCoords=True)
        if self.verbose_l2:
            print(f"{self.compound_code} NumAtomswithHs: {m.GetNumAtoms()}\n")
        # Embed the 3D structure

        # Write molecule
        if output_dir:
            mol_path = os.path.join(output_dir, f"{self.id}_Hs.mol")
            self.prepared_ligand_path = mol_path
            self.prepared_ligand_filename = os.path.basename(mol_path)
            out = Chem.SDWriter(mol_path)
            out.write(m)
            out.close()

        return m

    def draw_highlighted_mol(
        self,
        atom_indices,
        filename=None,
        img_size=(400, 400),
        save_3d_img=True
        ):
        """
        Draw the pose's rdkit.Chem.Mol with highlighted atoms.
        """
        filename3d = filename.replace(".png", "_3d.png")
        bond_indices = []
        for i in range(len(atom_indices)):
            for j in range(i + 1, len(atom_indices)):
                bond = self.ligand_mol.GetBondBetweenAtoms(
                    atom_indices[i], atom_indices[j]
                )
                if bond is not None:
                    bond_indices.append(bond.GetIdx())

        ligand_mol2d = Chem.Mol(self.ligand_mol)  # Copy the molecule
        ligand_mol2d.Compute2DCoords()

        img = Draw.MolToImage(
            ligand_mol2d,
            highlightAtoms=atom_indices,
            highlightBonds=bond_indices,
            size=img_size,
        )

        if filename:
            img.save(filename)
            if save_3d_img:
                img3d = Draw.MolToImage(
                    self.ligand_mol,
                    highlightAtoms=atom_indices,
                    highlightBonds=bond_indices,
                    size=img_size,
                )
                img3d.save(filename3d)
        else:
            img.show()

	# FIXME: Don't call self.mol here - find alternative
    def draw_mol3d(self, filename=None):
        """
        Draw the pose's rdkit.Chem.Mol. This is drawn considering 3D coords.
        """
        if filename:
            Draw.MolToFile(self.mol, filename)
        else:
            img = Draw.MolToImage(self.mol)
            img.show()

    def draw_mol(self, filename=None):
        """
        Draw the pose's rdkit.Chem.Mol. This is drawn considering 2D coords.
        """
        if filename:
            m = MolFromSmiles(self.smiles)
            Draw.MolToFile(m, filename)
        else:
            m = MolFromSmiles(self.smiles)
            img = Draw.MolToImage(m, size=(300, 300))
            img.show()

    def draw(self, inspirations=True, protein=False, **kwargs):
        """Render this pose (and its inspirations)"""

        if protein:
            from molparse.py3d import render

            return render(self.complex_system, **kwargs)

        mols = [self.mol]
        if inspirations and self.inspirations:
            mols += [i.mol for i in self.inspirations]

        return draw_mols(mols)

    @property
    def duck_feature_names(self):
        """Return the pose's protein feature names in a list format suitable for the DUck output"""
        # Combine only the second index of the dict valies into a LIST for DUck
        return [v[1] for v in self.fingerprint.values()]
        # return self.fingerprint_ext

    @property
    def protein_features(self):
        """Return the pose's protein features as a list"""
        return [v[0] for v in self.fingerprint_ext.values()]

    @property
    def ligand_features(self):
        """Return the pose's ligand features as a list"""
        return [v[1] for v in self.fingerprint_ext.values()]

    @property
    def mol(self):  # TODO: compare this to the mol property in pose.py
        """Returns a pose's rdkit.Chem.Mol"""
        if not self._mol:
            if self.mol_path.endswith(".mol"):
                if self.verbose:
                    print("")
                    logger.reading(self.mol_path)
                mol = mp.parse(self.mol_path, verbosity=False)

                if not mol:
                    logger.error(
                        f"[{self}] Error computing RDKit Mol from .mol={self.path}"
                    )
                    raise InvalidMolError

                self.mol = mol

            else:
                raise NotImplementedError

            if not mol:
                logger.error(f"Could not parse {self}.path={self.mol_path}")

        return self._mol

    @mol.setter
    def mol(self, m):
        """Set the pose's rdkit.Chem.Mol"""
        assert m
        self._mol = sanitise_mol(m)
        # self.db.update(table='pose', id=self.id, key='pose_mol', value=m.ToBinary())

    @property
    def features(self):
        """Returns the pose's features"""
        return mp.rdkit.features_from_mol(self.mol)


class InvalidMolError(Exception): ...


# Test
if __name__ == "__main__":
    from FragFeatures.fragalysis.target_parser import TargetParser

    target = TargetParser(
        "/Users/nfo24278/Documents/dphil/diamond/DuCK/structures/CHIKV_Mac"
    )
    lig = Pose(target.target_dir, "cx0270a")
    # lig.calculate_fingerprint()
    # print(lig.duck_feature_names)
    # print(target.target_dir)
    # print(lig.protein_path)
    # print(lig.mol_path)
    print("\n")
    lig.calculate_prolif_fp(
        "/Users/nfo24278/Documents/dphil/diamond/DuCK/code/features/Protein_preparation",
        rdkit_protein=True,
        verbose=True,
    )
