import freesasa as fs

pdb_file="/home/joaov/github/mmpbsa/testing/teste_gro.pdb" #frame00100.pdb"

fs.setVerbosity(1)
print(fs.getVerbosity())

structure = fs.Structure(pdb_file)
result =fs.calc(structure,fs.Parameters({'algorithm' : fs.ShrakeRupley,
                                               'n-points' : 10000}))
#result =fs.calc(structure,fs.Parameters({'algorithm' : fs.LeeRichards,
#                                               'n-slices' : 100}))

#area_classes = fs.classifyResults(result, structure)
print("Total : %.2f A2" % result.totalArea())
#for key in area_classes:
#    print(key, ": %.2f A2" % area_classes[key])
#    

selections = fs.selectArea(('all, resi 1-198','chainA, resi 1-99', 'chainB, resi 100-198'), structure, result)
list_keys=[]
for key in selections:
    print (key, ": %.2f A2" % selections[key])
    list_keys.append(selections[key])
    # Saving the objects:
g = 0.00542 # gamma -> kcal/(mol‚A2)
b = 0.92    # beta  -> kcal/mol

P=list_keys[0]*g+b
A=list_keys[1]*g+b
B=list_keys[2]*g+b

energy= P - (A + B)
print(P,A,B)
print("energy=",energy*4.184)
