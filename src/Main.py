from Structure import Structure
#from ConformerGenerator import ConformerGenerator
#from Analyzer import Analyzer



def main():
    # example workflow
    #Molecule = Structure()
    #Molecule.get_structure("../.xyz-Inputdateien/Tryptophan.xyz")
    #Generator = ConformerGenerator()
    #Generator.generate_conformers(Molecule)

    Molecule = Structure()
    Molecule.get_structure("/home/dario/SUPRA-conformer/inputfiles/Tyrosin.xyz")
    Molecule.create_mol_file()



if __name__ == "__main__":
    main()
