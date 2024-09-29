from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import pandas as pd
import os


#Load CSV file containing compounds of interest, requires that at least one column contains the compound name, and another contains the SMILES
df = pd.read_csv("compounds.csv") 
smiles_column = "Smiles" #Adjust as needed
id_column = "Name"  #Adjust as needed

#Prepare directory for batch processing of compounds separating by every 1000 compounds
output_base_dir = "sdf_files"
compounds_per_folder = 1000  #Set limit for compounds per folder
folder_count = 1
compound_count = 0

#Generate the first folder
current_output_dir = os.path.join(output_base_dir, f"batch_{folder_count}")
if not os.path.exists(current_output_dir):
    os.makedirs(current_output_dir)

for index, row in df.iterrows():
    smiles = row[smiles_column]
    compound_name = str(row[id_column]) 

    #Convert a molecules Smiles into a mol object in RdKit and add hydrogens
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    #Generate 2D coordinates of the mol object
    AllChem.Compute2DCoords(mol)

    #Generate 3D conformation of the mol object and optimize using the MMF94S force-field
    AllChem.EmbedMolecule(mol, AllChem.ETKDG()) 
    AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s")  

    #Save the mol obj as a .sdf into the output directory
    sdf_filename = os.path.join(current_output_dir, f"{compound_name}.sdf")
    writer = Chem.SDWriter(sdf_filename)
    writer.write(mol)
    writer.close()

    #Increment through the compound list
    compound_count += 1

    # After every n compounds, create a new folder
    if compound_count % compounds_per_folder == 0:
        folder_count += 1
        current_output_dir = os.path.join(output_base_dir, f"batch_{folder_count}")
        if not os.path.exists(current_output_dir):
            os.makedirs(current_output_dir)

print("Compound Conversion Complete")



