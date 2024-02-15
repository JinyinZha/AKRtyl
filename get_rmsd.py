import os
import prody
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def find_conf(plot,cv1,cv2):   
    confid = []
    x = []
    y = []
    points = np.array([cv1,cv2]).T
    for i in range(1,len(plot[1])-2):
        for j in range(1,len(plot[2])-2):
            if plot[0][i,j] == 0:
                continue
            if plot[0][i-1,j-1] == 0 or plot[0][i-1,j] == 0 or plot[0][i-1,j+1] == 0 or plot[0][i,j-1] == 0  or plot[0][i,j+1] == 0 or plot[0][i+1,j-1] == 0 or plot[0][i+1,j] == 0 or plot[0][i+1,j+1] == 0:
                continue
            if plot[0][i,j] >= plot[0][i-1,j-1] and plot[0][i,j] >= plot[0][i-1,j] and plot[0][i,j] >= plot[0][i-1,j+1] and plot[0][i,j] >= plot[0][i,j-1] and plot[0][i,j] >= plot[0][i,j+1] and plot[0][i,j] >= plot[0][i+1,j-1] and plot[0][i,j] >= plot[0][i+1,j] and plot[0][i,j] >= plot[0][i+1,j+1]:
                st1 = plot[1][i]
                ed1 = plot[1][i+1]
                st2 = plot[2][j]
                ed2 = plot[2][j+1]
                d2center = np.linalg.norm(points - np.array([0.5*(st1+ed1),0.5*(st2+ed2)]),axis=1)
                min_id = np.argwhere(d2center == np.min(d2center)).flatten()[0]                
                confid.append(min_id)
                x.append(cv1[min_id])
                y.append(cv2[min_id])
    return confid,x,y

#1, Load Data
name = os.getcwd().split("/")[-1] 
pdb = prody.parsePDB("%s_md.pdb"%(name))
pdb0 = pdb.copy()
rmsds = []
chains = "A"
rs = [list(range(0,50000)),
      list(range(0,50000)),
      list(range(0,50000)),
]
for i in range(1,4,1):
    traj = prody.parseDCD("%s_%d.dcd"%(name,i))
    traj.setAtoms(pdb0)
    pdb.setCoords(traj.getCoordsets())
    r = rs[i-1]
    ca_coord0 = prody.parsePDB("/home/jyzha/project/AKR/ndp_mono/ndp_mono_md.pdb").select("name CA").getCoords()
    for c in chains:
        e = prody.Ensemble()
        e.setCoords(ca_coord0)
        e.addCoordset(pdb.select("name CA and chain %s"%(c)).getCoordsets())
        e.superpose()
        rmsds.append(e.getRMSDs())
     
np.save("%s_rmsd.npy"%(name),rmsds)
