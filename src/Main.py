from Structure import Structure
from revConformerGenerator import ConformerGenerator
#from ConformerGenerator import ConformerGenerator
#from Analyzer import Analyzer



def main():
    # example workflow
    #Molecule = Structure()
    #Molecule.get_structure("../.xyz-Inputdateien/Tryptophan.xyz")
    #Generator = ConformerGenerator()
    #Generator.generate_conformers(Molecule)
    
    filename = "home/baum/SUPRA-conformer/inputfiles/Alanin.xyz"
    molecule = Structure(filename)
    



if __name__ == "__main__":
    main()
