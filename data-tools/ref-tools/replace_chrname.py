import sys
fin = open(sys.argv[1])
fout = open(sys.argv[1] + ".diffcontigname", "w")
data = fin.readlines()
for line in data:
    if line[0] == '>':
        tokens = line.strip().split(" ")
        if tokens[0] == ">gi|224589800|ref|NC_000001.10|":
            fout.write(">1 " + " ".join(tokens[1:]) + "\n")            
        if tokens[0] == ">gi|224589811|ref|NC_000002.11|":
            fout.write(">2 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589815|ref|NC_000003.11|":
            fout.write(">3 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589816|ref|NC_000004.11|":
            fout.write(">4 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589817|ref|NC_000005.9|":
            fout.write(">5 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589818|ref|NC_000006.11|":
            fout.write(">6 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589819|ref|NC_000007.13|":
            fout.write(">7 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589820|ref|NC_000008.10|":
            fout.write(">8 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589821|ref|NC_000009.11|":
            fout.write(">9 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589801|ref|NC_000010.10|":
            fout.write(">10 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589802|ref|NC_000011.9|":
            fout.write(">11 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589803|ref|NC_000012.11|":
            fout.write(">12 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589804|ref|NC_000013.10|":
            fout.write(">13 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589805|ref|NC_000014.8|":
            fout.write(">14 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589806|ref|NC_000015.9|":
            fout.write(">15 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589807|ref|NC_000016.9|":
            fout.write(">16 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589808|ref|NC_000017.10|":
            fout.write(">17 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589809|ref|NC_000018.9|":
            fout.write(">18 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589810|ref|NC_000019.9|":
            fout.write(">19 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589812|ref|NC_000020.10|":
            fout.write(">20 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589813|ref|NC_000021.8|":
            fout.write(">21 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589814|ref|NC_000022.10|":
            fout.write(">22 " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589822|ref|NC_000023.10|":
            fout.write(">X " + " ".join(tokens[1:]) + "\n")
        if tokens[0] == ">gi|224589823|ref|NC_000024.9|":
            fout.write(">Y " + " ".join(tokens[1:]) + "\n")
    else:
        fout.write(line)
fin.close()
fout.close()
