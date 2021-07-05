import sys

if len(sys.argv) != 2:
    print("usage: ", sys.argv[0], " file ")
    exit(1)

names = ["DE_int", "DE_tildePauli", "DEexc", "DEpauli", \
        "DEpauli", "DEelect", "DEorb"]
energies = {}
for name in names:
    energies[name] = 0.0


fp = open(sys.argv[1])

for line in fp:

    for name in names:
        if line.find(name) >= 0:
            energies[name] = float(line.split()[-1])
            
    for i in range(1, 13, 2):
        if line.find("nocv_%d pair eigenvalue"%(i)) >= 0:
            energies["nocv_%d"%(i)] = 2.0 * float(line.split()[-1])
            

print("\\begin{table}[]")
print("\\begin{tabular}{ll}")
for name in energies:
    print(name.replace("_", "") , " & ", energies[name] , "\\\\")
print("\\end{tabular}")
print("\\end{table}")
print("\\end{document}")


