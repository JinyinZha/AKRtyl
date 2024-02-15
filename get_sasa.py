from pymol import cmd
import prody
import numpy as np
import copy
import os

name = os.getcwd().split("/")[-1]
pdb = prody.parsePDB("%s_md.pdb"%(name))

chains = "A"
cmd.delete("all")
cmd.set("dot_solvent",1)
cmd.set("dot_density",4)


surface_id_li = [21, 22, 25, 26, 50, 51, 52, 53, 55, 56, 59, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 101, 102, 130, 131, 132, 160, 161, 188, 214, 215, 218, 219, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244]
poc_sels = [" or ".join(["resi %d"%(i) for i in pdb.select("chain %s and name CA"%(c)).getResnums()[surface_id_li]]) for c in chains]
print(poc_sels)


all_sasas = []
for n in range(1,4,1):
    traj = prody.parseDCD("%s_%d.dcd"%(name,n))
    traj.setAtoms(pdb)
    print("traj %d"%(n))
    sasas = []
    for i in range(len(chains)):
        this_sasas = []
        poc_sel = poc_sels[i]
        '''
        cmd.delete("all")
        cmd.load("../ndp_sub_mono/ndp_sub_mono_md.pdb")
        cmd.load("%s_md.pdb"%(name))
        cmd.remove("%s_md and not chain %s"%(name,chains[i]))
        cmd.remove("%s_md and not polymer"%(name))
        cmd.align(mobile="ndp_sub_mono_md",target="%s_md"%(name),cycles=50)
        cmd.remove("ndp_sub_mono_md and not resname TYK")
        cmd.save("tmp.pdb","byres resname TYK around 6")
        tmp = prody.parsePDB("tmp.pdb")
        poc_sel_li = []
        for ri in np.unique(tmp.getResnums()):
            poc_sel_li.append("resi %d"%(ri))
        poc_sel = " or ".join(poc_sel_li)
        print("Chain %s: sel=%s"%(chains[i],poc_sel))
        '''
        
        for n in range(0,50000,100):
            frame = traj[n]
            prody.writePDB("tmp.pdb",frame)
            cmd.delete("all")
            cmd.load("tmp.pdb")
            cmd.remove("not chain %s"%(chains[i]))
            cmd.remove("not polymer")
            sasa = cmd.get_area(poc_sel)
            this_sasas.append(sasa)
            print(n,chains[i],n,sasa,np.mean(this_sasas))
        sasas.append(copy.deepcopy(this_sasas))
        print("----------------------------------------------------")
    all_sasas.append(copy.deepcopy(sasas))
    print("----------------------------------------------------")
    print("----------------------------------------------------")
np.save("%s_pocSASA.npy"%(name),all_sasas)

