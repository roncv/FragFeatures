"""
Protein preparation module.
"""

import logging

# import parmed
import os
import shutil

from openmm import VerletIntegrator
from openmm.app import PME, ForceField, Modeller, PDBFile, Simulation
from openmm.unit import picoseconds

from FragFeatures.utils import timefunction

logger = logging.getLogger("FragFeatures")


class ProteinPreparation:
    """
    Prepare a protein for simuilation - add hydrogens.
    """

    def __init__(
        self,
        protein_path,
        output_dir,
        protein_id=None,
        minimize=True,
        pH=7.8,
        verbose=False,
        verbose_l2=False,
    ):
        self.protein_path = protein_path
        self.protein_file = os.path.basename(self.protein_path)
        self.output_dir = output_dir
        self.protein_id = protein_id
        self.protein_output_path = os.path.join(
            self.output_dir, os.path.basename(self.protein_path)
        )
        self.minimize = minimize
        self.pH = pH
        self.termini = True
        self.verbose = verbose
        self.verbose_l2 = verbose_l2

    def create_output_dir(self):
        """
        Create the output directory for the prepared protein.
        """
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)
        # Duplicate the protein file
        shutil.copy(self.protein_path, self.protein_output_path)

    @timefunction
    def prepare_protein(self):
        """
        Prepare the protein and place it in the output directory.
        """
        self.create_output_dir()
        self.add_hydrogens(minimize=self.minimize, pH=self.pH)

    def add_hydrogens(self, minimize=True, pH=7.8):
        """
        Add hydrogens to the protein.
        """
        # logger.info(f'Running command: {command}')
        residues_to_remove = [
            "DMS",
            "TRS",
            "LIG",
            "CL",
        ]  # TODO: Make into constant # +'CL'
        self.prepared_protein_path = self.protein_output_path.replace(
            ".pdb", "_prepared.pdb"
        )
        self.protein_filename = os.path.basename(self.prepared_protein_path)

        shutil.copy(self.protein_output_path, self.prepared_protein_path)

        pdb = PDBFile(self.protein_output_path)
        forcefield = ForceField("amber99sb.xml", "tip3p.xml")
        modeller = Modeller(pdb.topology, pdb.positions)

        for resname in residues_to_remove:
            modeller.delete(
                [
                    residue
                    for residue in modeller.topology.residues()
                    if residue.name == resname
                ]
            )
        modeller.deleteWater()

        openmm_resids = [residue.index for residue in modeller.topology.residues()]
        # Debugging
        if self.verbose_l2:
            unique_residues = set(
                [residue.name for residue in modeller.topology.residues()]
            )
            print(f"Unique Residues: \n\n{unique_residues}")
            resnames = [residue.name for residue in modeller.topology.residues()]
            resids = [residue.id for residue in modeller.topology.residues()]
            residues = [
                f"{resname}{resid}_{residx+1}"
                for resname, resid, residx in zip(resnames, resids, openmm_resids)
            ]
            print(f"\nResidues: \n\n{residues}")

        if self.verbose:
            if self.protein_id:
                print(f"\nPreparing {self.protein_file} protein...")
            print("\nAdding hydrogens...")

        try:
            modeller.addHydrogens(forcefield, pH=pH)
        except:
            # Try removing the termini
            print("Fixing protein...")
            print("Removing termini...")
            if self.verbose:
                print(
                    f"Termini: (lower_idx: {min(openmm_resids)}, higher_idx: {max(openmm_resids)})"
                )
            if self.verbose_l2:
                print(f"\nResIdxs:\n{openmm_resids}\n")
            termini_idxs = [min(openmm_resids), max(openmm_resids) - 1]
            for residx in termini_idxs:
                modeller.delete(
                    [
                        residue
                        for residue in modeller.topology.residues()
                        if residue.index == residx
                    ]
                )
            openmm_resids_cut = [
                residue.index for residue in modeller.topology.residues()
            ]
            if self.verbose_l2:
                print(f"\nResIdxs Cut:\n{openmm_resids_cut}\n")
            modeller.addHydrogens(forcefield, pH=pH)
            self.termini = False

        system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)
        integrator = VerletIntegrator(0.001 * picoseconds)
        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)

        if minimize:
            if self.verbose:
                print("Minimizing...")
            simulation.minimizeEnergy(maxIterations=100)

        if self.verbose:
            print(f"Saving {self.protein_filename}...")
        positions = simulation.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(
            simulation.topology, positions, open(self.prepared_protein_path, "w")
        )

        if self.verbose:
            print("Done.\n")

    def get_prepared_protein_path(self):
        """
        Return the path to the prepared protein.
        """
        return self.prepared_protein_path


    def get_protein_filename(self):
        """
        Return the name of the protein file.
        """
        return self.protein_filename

    def cleanup(self):
        """
        Remove the prepared protein.
        """
        os.remove(self.prepared_protein_path)


# Test
if __name__ == "__main__":
    protein = ProteinPreparation(
        protein_path="/Users/nfo24278/Documents/dphil/diamond/DuCK/structures/CHIKV_Mac/aligned_files/cx0270a/cx0270a_apo.pdb",
        output_dir="/Users/nfo24278/Documents/dphil/diamond/DuCK/code/features/Protein_preparation",
    )
    protein.prepare_protein()
    print(protein.get_prepared_protein_path())
