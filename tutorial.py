import bpRNAStructure.StructureCreator as STC

# call class method to create new structure from file: bpRNA_CRW_2.st
structure = STC.StructureCreator.create("./bpRNA_CRW_2.st")

print(structure)  # RNA: bpRNA_CRW_2

# iterate through all stems in the molecule and print their sequence
for stem in structure.stems():
    print(stem.sequence())
