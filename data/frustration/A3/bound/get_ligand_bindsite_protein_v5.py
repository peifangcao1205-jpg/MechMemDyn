import numpy as np
import sys
import math
import os
residue_type = ["ARG","LYS","ASN","GLN","GLU","ASP","HIS","TYR","TRP","SER","THR","GLY","PRO","ALA","MET","CYS","PHE","LEU","VAL","ILE"]
def vector(p1, p2):
    return [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]

def vabs(a):
    return math.sqrt(pow(a[0],2)+pow(a[1],2)+pow(a[2],2))

def angle(p1,p2):
    return math.acos((p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2])/vabs(p1)/vabs(p2))
def get_protein_ligand_allatom(pdbfile):
    id_res_old = -999
    atom = []
    res = []
    protein = []
    chain_id = []
    res_id = []
    ligand_atom = []
    ligand = []
    ligand_name = []
    name_ligand_old = "HOH"
    with open(pdbfile,"r") as fopen:
         lines = fopen.readlines()
    for line in lines:
        if len(line.split()) > 5:
           if line[0:4] == "ATOM":
              id_res = int(line[22:26])
              #chain_id = line[21:22]
              #print id_res_old
              if id_res != id_res_old:
                 if id_res_old != -999:
                    protein.append(res)
                    chain_id.append(line[21:22])
                    res_id.append(int(line[22:26]))
                    #print id_res_old
                 id_res_old = id_res
                 res = []
                 x=float(line[30:38])
                 y=float(line[38:46])
                 z=float(line[46:54])
                 atom = [x,y,z]   
                 res.append(atom)
                 #if line[13:17] == "CA":
                 #   chain_id.append(line[21:22])
                 #   res_id.append(int(line[22:26]))
              else:
                 x=float(line[30:38])
                 y=float(line[38:46])
                 z=float(line[46:54])
                 atom = [x,y,z]
                 res.append(atom)
                 #if line[13:17] == "CA":
                 #   chain_id.append(line[21:22])
                 #   res_id.append(int(line[22:26]))
           if line[0:6] == "HETATM":
                 name_ligand = line[17:20]
                 if name_ligand != 'HOH':
                       x=float(line[30:38])
                       y=float(line[38:46])
                       z=float(line[46:54])
                       ligand_atom.append([x,y,z])
    protein.append(res)
    #print ligand_atom
    #print chain_id
    return ligand_atom,protein,chain_id,res_id

def compute_bindsite(ligand_atom,protein,chain_id,res_id):
    data = ""
    cutoff = 6
    dist = 999
    #print len(protein)
    for i in range(len(protein)):
        for j in range(len(protein[i])):
            for k in range(len(ligand_atom)):
                  dist = vabs(vector(protein[i][j],ligand_atom[k]))
                  if dist < cutoff:
                      break
            if dist < cutoff:
                data += chain_id[i] + str(res_id[i])+" "
                break
    return data 

def main():
    pdbfile = sys.argv[1]
    ligand_atom,protein,chain_id,res_id = get_protein_ligand_allatom(pdbfile)
    data = compute_bindsite(ligand_atom,protein,chain_id,res_id)
    #name = pdbfile.split("/")[1].split(".")[0] + ".list"
    name = "bres.list"
    with open(name,"w") as fwrite:
         fwrite.writelines(data) 

if __name__ == '__main__':
    main()
  
        

                        
               
    
