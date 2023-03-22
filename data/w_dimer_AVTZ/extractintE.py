import sys

if len(sys.argv) != 2:
    print("usage: ", sys.argv[0], " file ")
    exit(1)

names = ["DE_int", "DE_tildePauli", "DEexc", "DEpauli", \
        "DEpauli", "DEelect", "DEorb"]
energies = {}
nocv     ={}


fp = open(sys.argv[1])

for line in fp:

    if line.find("DE_int") >= 0:
            aa = "$\Delta E_{int}$"
            energies.update({aa : float(line.split()[-1])})

    if line.find("DE_tildePauli") >= 0:
            aa = "$\Delta E_{ETS-Pauli}$"
            energies.update({aa : float(line.split()[-1])})

    if line.find("DEexc") >= 0:
            aa = "$\Delta E_{xc}$"
            energies.update({aa : float(line.split()[-1])})

    if line.find("DEpauli") >= 0:
            aa = "$\Delta E_{Pauli}$"
            energies.update({aa : float(line.split()[-1])})

    if line.find("DEelect") >= 0:
            aa = "$\Delta E_{elstat}$"
            energies.update({aa : float(line.split()[-1])})

    if line.find("DEorb") >= 0:
            aa = "$\Delta E_{orb}$"
            energies.update({aa : float(line.split()[-1])})
    
    j = 1

    for i in range(1, 13, 2):
        if line.find("nocv_%d pair eigenvalue"%(i)) >= 0:
            aa = "nocv_%d"%(i)
            aa = "$\Delta E^{%d}_{orb}$"%(i)
            nocv.update({aa : [float(line.split()[6]), float(line.split()[3])]})
            j += 1




print("\\begin{table}[]")
print("\\caption{}")
print("\\label{}")
print("\\begin{tabular}{ll}")
print("\\hline")
for name in energies:
    print(name.ljust(25) + " & %12.2f"%( energies[name]*627.51 ) + "   \\\\")

for name in nocv:
    print(name.ljust(25) + " & %12.2f"%( nocv[name][0]*627.51*2 )+"(%12.4f)"%(abs(nocv[name][1])*2) + "   \\\\")
print("\\end{tabular}")
print("\\end{table}")
#print("\\end{document}")

