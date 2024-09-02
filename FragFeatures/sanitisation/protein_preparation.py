"""
Protein preparation module.
"""

from FragFeatures.utils import timefunction, dict_to_json

import os
import json
import shutil
import subprocess

import logging
logger = logging.getLogger('FragFeatures')



class ProteinPreparation:
    """
    Prepare a protein for simuilation - add hydrogens.
    """
    def __init__(self, protein_path, output_dir):
        self.protein_path = protein_path
        self.output_dir = output_dir
        self.protein_output_path = os.path.join(self.output_dir, os.path.basename(self.protein_path))

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
        self.add_hydrogens(minimize=True)


    def add_hydrogens(self, minimize=True, pH=7.8):
        """
        Add hydrogens to the protein.
        """
        # logger.info(f'Running command: {command}')
        residues_to_remove = ['DMS', 'TRS', 'LIG', 'CL'] # TODO: Make into constant
        self.prepared_protein_path = self.protein_output_path.replace('.pdb', '_prepared.pdb')

        shutil.copy(self.protein_output_path, self.prepared_protein_path)

        # Using OpenMM to add hydrogens
        from openmm.app import PDBFile, Modeller, ForceField, Simulation, PME
        from openmm import VerletIntegrator
        from openmm.unit import picoseconds

        pdb = PDBFile(self.protein_output_path)
        forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
        modeller = Modeller(pdb.topology, pdb.positions)

        # Debugging
        all_residues = set([residue.name for residue in modeller.topology.residues()])
        print(all_residues)

        for resname in residues_to_remove:
            modeller.delete([residue for residue in modeller.topology.residues() if residue.name == resname])
        modeller.deleteWater()

        print('Adding hydrogens...')
        modeller.addHydrogens(forcefield, pH=pH)
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)
        integrator = VerletIntegrator(0.001*picoseconds)
        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)

        if minimize:
            print('Minimizing...')
            simulation.minimizeEnergy(maxIterations=100)

        print('Saving...')
        positions = simulation.context.getState(getPositions=True).getPositions()
        PDBFile.writeFile(simulation.topology, positions, open(self.prepared_protein_path, 'w'))
        print('Done')



    def get_prepared_protein_path(self):
        """
        Return the path to the prepared protein.
        """
        return self.prepared_protein_path


    def cleanup(self):
        """
        Remove the prepared protein.
        """
        os.remove(self.prepared_protein_path)


# Test
if __name__ == '__main__':
    protein = ProteinPreparation(protein_path='/Users/nfo24278/Documents/dphil/diamond/DuCK/structures/CHIKV_Mac/aligned_files/cx0270a/cx0270a_apo.pdb',
                                 output_dir='/Users/nfo24278/Documents/dphil/diamond/DuCK/code/features/Protein_preparation'
                                 )
    protein.prepare_protein()
    print(protein.get_prepared_protein_path())
