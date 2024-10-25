from Bio.PDB import *
from scipy.spatial import distance_matrix
from scipy.spatial.transform import Rotation as Rot

import os
import numpy as np
import pandas as pd
import copy

BOND_CUT = 2
INTERACT = 4.2

a1 = 1.5
a2 = 8.5

VDW_CUT = 1.75

INTERACT_RING = 5

#Determine which ones are glycosylated and which are non-covalent

def get_ligand_coor(structure):

    coor_c = []
    coor_p = []
    res_c = []
    res_p = []

    coor_x = []
    res_x = []

    models = structure.get_models()
    models = list(models)
    for m in range(len(models)):
        chains = list(models[m].get_chains())
        for c in range(len(chains)):
            residues = list(chains[c].get_residues())
            for r in range(len(residues)):
                res = residues[r].get_resname()
                if res == 'HOH':
                    continue;

                atoms = list(residues[r].get_atoms())

                for a in range(len(atoms)):
                    at = atoms[a]

                    if 'H' in at.get_name():
                        continue;

                    #print(str(residues[r].get_parent().id).strip())

                    #print(str(residues[r].get_resname()))
                    if 'LG1' in str(residues[r].get_resname()) or 'LIG' in str(residues[r].get_resname()): #or 'UNK' in str(residues[r].get_resname()):
                        coor_c.append( at.get_coord() )
                        res_c.append( [ str(residues[r].id[1]).strip(), str(chains[c].id).strip(), str(residues[r].get_resname()), str(at.get_name()) ] )
                        #res_c.append( [ str(residues[r].id[1]).strip(), str(residues[r].get_resname()) ] )
                    else:
                        coor_p.append( at.get_coord() )
                        res_p.append( [ str(residues[r].id[1]).strip(), str(chains[c].id).strip(), str(residues[r].get_resname()), str(at.get_name()) ] )
                        #res_p.append( [ str(residues[r].id[1]).strip(), str(residues[r].get_resname()) ] )
                    #else:
                    #    coor_x.append( at.get_coord() )
                    #    res_x.append( [ str(residues[r].id[1]).strip(), str(chains[c].id).strip(), str(residues[r].get_resname()) ] )
                        #res_x.append( [ str(residues[r].id[1]).strip(), str(residues[r].get_resname()) ] )


    return res_p, coor_p, res_c, coor_c

class glycan():
    """
    Class object for a GLYCAN

    Args:
        coors (arr nx3): coordinates of heavy atoms
        atom_names (arr str): names of the atoms

    Variables:
        name, coor, atom_names
        adj_mat (nxn): One-hot of bonded atoms
        edges (nx?): array of arrays of the non-sparse edge connections
        ring_atom (arr nx[5,6]x1 or nx6x1): defines which atoms are in the ring
    """

    def __init__(self,coor,atom_names,BOND_CUTOFF=1.75):

        self.coor = coor
        #print(self.coor)
        self.atom_names = atom_names

        self.BOND_CUTOFF = BOND_CUTOFF


        #initialize empty variables
        self.adj_mat = []
        self.edges = []
        self.ring_atom = []
        self.com = []
        self.ring_atom_plus = []

        self.calc_adjacency()

        ope = []

        for jj in range(len(self.coor)):
            o = self.calc_ring(jj)
            ope.append(o)
            self.calc_adjacency()

        ring = []
        ring.append(ope[0])
        for jj in range(1,len(ope)):
            if type(ope[jj]) == bool:
                continue;

            skip = False
            for kk in range(len(ring)):
                if ring[kk][0] == ope[jj][0]:
                    skip = True
            if skip:
                continue;

            ring.append(ope[jj])


        ring_plus = []
        for jj in range(len(ring)):
            ring_plus.append([])
            r = []
            for kk in ring[jj]:
                r.append(self.coor[kk])
            #print(self.coor[kk])
            #print(self.coor)
            d = distance_matrix(r,np.array(self.coor))
            d = d < self.BOND_CUTOFF
            #print(np.shape(d))
            d = np.sum(d,axis=0) >= 1
            #print(np.shape(d))
            for ll in range(len(d)):
                if d[ll] and ll not in ring_plus[jj]:
                    ring_plus[jj].append(ll)
        self.ring_atom_plus = ring_plus

        #for jj in range(len(ring)):
        #    print('a',ring[jj],'\n','p',ring_plus[jj])


        self.ring_atom = ring
        self.ring_atom_name, self.ring_com = self.get_ring_atom_name()
        #self.print_variables()

    def calc_adjacency(self):
        #get the adjacency matrix and edge list of the carb

        #calculate atom-atom distances and set cutoffs
        dm = distance_matrix(self.coor,self.coor)
        #print(dm)
        adj_mat = dm < self.BOND_CUTOFF;
        #no self interactions
        for i in range(len(adj_mat)):
            adj_mat[i,i] = 0

        #get the list of the adjacency matrix
        edge_list = [];
        for ii in range(len(adj_mat)):
            edge_list.append([])
            for jj in range(len(adj_mat)):
                if adj_mat[ii,jj]:
                    edge_list[ii].append(jj)

        #store local variables into class variables
        self.adj_mat = adj_mat
        self.edges = edge_list
        return

    #recursive algo to get cycle of the graph
    def visit(self,n,edge_list,visited,st):
        """
        Args:
            n - node we are searching from
            edge_list - adjacency of each node, is periodically
                modified to remove connection to parent coming from
            st - start node
        Returns:
            arr - array of the cycle found
        """
        #print(edge_list)
        #print(n)

        if n == st and visited[st] == True:
            return [n]

        visited[n] = True
        r = False
        arr = []
        #print(n,edge_list[n],visited)

        for e in edge_list[n]:
            #if n in edge_list[e]:
            try:
                edge_list[e].remove(n)
            except:
                continue;
            #print('\t',e)

            r = self.visit(e,edge_list,visited,st)
            #print('\t\t',r)

            if type(r) != bool:
                arr.append(n)
                for j in r:
                    arr.append(j)

        if arr == []:
            return False

        return arr

    def calc_ring(self,i):
        #gets the ring atoms, calls recursive visit function
        ring = self.visit(i,copy.deepcopy(self.edges),np.zeros(len(self.coor)),i)
        ind = 0
        while type(ring) == bool:
            #print('ringboi',ind,len(self.coor))
            ring = self.visit(ind,copy.deepcopy(self.edges),np.zeros(len(self.coor)),ind)
            ind += 1;
            if ind >= len(self.coor):
                break;

        #print(ring)
        self.ring_atom = np.unique(ring).astype(int)

        return self.ring_atom

    def get_ring_atom_name(self):
        #gets the ring_atom_names in PDB notation and the com of each ring
        r = []
        com = []
        for jj in self.ring_atom:
            r.append([])
            com.append(np.array([0.,0.,0.]))
            for kk in jj:
                r[-1].append(self.atom_names[kk])
                com[-1] += np.array(self.coor[kk])
            com[-1] /= len(r[-1])
        return r, np.array(com)

def get_rings(file):
    """
    input:
        file (str): file name string
    output:
        ring_atom_name (arr, str): PDB names of glycan ring ATOMS
        ring_com (arr, float): Center of Mass (COM) of all ring atoms
        gly (class): raw glycan class for further analysis if needed
    """
    parser=PDBParser()
    structure=parser.get_structure("prot", file)
    _,_, res, coor = get_ligand_coor(structure)
    at = []
    for ii in res:
        at.append(ii[-1])
    #print(at)
    gly = glycan(coor,at)
    return gly.ring_atom_name, gly.ring_com, gly

def find_interactChains(coor_c,coor_p,res_c,res_p,INTERACT=4.2):
    #determine chain-chain interactions
    d = distance_matrix(coor_c,coor_p) < INTERACT
    a = np.array( np.where(d == 1) )
    a = np.array(a)

    #chain_int = {}
    chain_int = []
    for ii in range(a.shape[1]):
        #res1 = res_c[ a[0,ii] ]
        res2 = res_p[ a[1,ii] ]

        #print(res1,res2)

        chain_int.append(res2)

    return chain_int

def find_interactRingRes(rcom,pc,pr,INTERACT=6.0):
    #determine chain-chain interactions
    d = distance_matrix(rcom,pc) < INTERACT

    cint = []
    for ii in range(len(rcom)):
        cint.append([])
        for jj in range(len(pr)):
            if d[ii,jj]:
                res2 = int(pr[jj][0])
                if res2 not in cint[-1]:
                    cint[-1].append(res2)
    return cint

def find_interactRingAtomRes(rcom,gly,pc,pr,INTERACT=6.0):
    #determine chain-chain interactions
    cint = []
    for jj in range(len(gly.ring_atom_plus)):

        #print(gly.ring_atom[jj])

        at = []
        for ii in gly.ring_atom_plus[jj]:
            at.append(gly.coor[ii])

        d = distance_matrix(at,pc) < INTERACT
        cint.append([])

        for ii in range(len(at)):

            for jj in range(len(pr)):
                if d[ii,jj]:
                    res2 = int(pr[jj][0])
                    if res2 not in cint[-1]:
                        cint[-1].append(res2)
    return cint

def fnat(cint,cint_):
    #gets Fnat - fraction of natural contacts
    f ,n = 0,0
    for ii in range(len(cint_)):
        for jj in range(len(cint_[ii])):
            n += 1
            if ii < len(cint):
                if cint_[ii][jj] in cint[ii]:
                    f += 1
    #print(f,n,f/n)
    return f / n

def adjusted_fnat(cint,cint_):
    #gets Fnat - fraction of natural contacts
    f = 0;
    for ii in range(-1,len(cint)):
        new_c = []
        for jj in range(len(cint)):
            new_c.append(cint[(jj + ii) % len(cint)])
        o = fnat(new_c,cint_)
        if o > f:
            f = o

    return f


def dice(wt_r,pred_r):
    y, y_hat = [], []

    for ii in wt_r:
        if ii[0] not in y:
            y.append(ii[0])
            #print(ii)
    for ii in pred_r:
        if ii[0] not in y_hat:
            y_hat.append(ii[0])
            #print('\t',ii)
    y = np.sort(np.array(y).astype(int))
    y_hat = np.sort(np.array(y_hat).astype(int))

    #print(y,'\n',y_hat)
    a = np.max(y_hat)
    if np.max(y) > a:
        a = np.max(y)

    y_arr = np.zeros(a + 500)
    y_pred_arr = np.zeros(a + 500)

    y_arr[y] = 1
    y_pred_arr[y_hat] = 1

    #print(y_arr * y_pred_arr)

    d = 2 * np.sum(y_arr * y_pred_arr) / (np.sum(y_arr) + np.sum(y_pred_arr))
    return d

def rms(x,y):
    if len(x) != len(y):
        print('lengths of ligands do not match')
        return -1
    r = 0
    for ii in range(len(x)):
        r += np.linalg.norm(x[ii] - y[ii]) ** 2
    r /= len(x)
    return np.sqrt(r)

def get_all_info(file):
    """
    input:
        file (str): file name string
    output:
        prot_res (arr, str): PDB names of protein CA atoms
        prot_coor (arr, float): coordinates of the CA atoms
        int_res (arr, str): PDB names of protein residues interacting with rings
        ring_atom_name (arr, str): PDB names of glycan ring ATOMS
        ring_com (arr, float): Center of Mass (COM) of all ring atoms
        gly_coor (arr,float): all atom coordinates of the glycan
        gly (class): raw glycan class for further analysis if needed
    """
    parser=PDBParser()
    structure=parser.get_structure("prot", file)
    pr,pc, res, coor = get_ligand_coor(structure)

    #get glycan_info
    at = []
    for ii in res:
        at.append(ii[-1])
    #print(at)
    gly = glycan(coor,at)

    #get protein info
    prot_res, prot_coor = [], []
    for ii in range(len(pr)):
        if pr[ii][-1] == 'CA':
            prot_res.append(pr[ii])
            prot_coor.append(pc[ii])

    #get interact info
    int = find_interactChains(coor,pc,res,pr)

    cint = find_interactRingAtomRes(gly.ring_com,gly,pc,pr,INTERACT = INTERACT_RING)

    return pr, pc, int, gly.ring_atom_name, gly.ring_com, gly.coor, gly, cint

def adjusted_lrms(gc,gc_,g,g_):

    r = 1e+10
    big_no = []

    for jj in range(len(gc_)):
        dm = distance_matrix(gc,gc_)
        no = []
        ope = []
        curr_no = []
        u_gc_ = []
        for kk in range(len(gc_)):
            ii = (jj + kk) % len(gc_)
            #print(ii)

            a = np.argmin(dm[:,ii])

            #print(dm[:,ii])

            skip = False
            curr_no = []
            if g.atom_names[a][0] == g_.atom_names[ii][0]:
                no.append(a)
                ope.append([a,ii])
                dm[:,ii] = 1e10
                dm[a,:] = 1e10
                skip = True
                u_gc_.append(gc_[ii])
            else:
                curr_no.append(a)

            cnt = 0
            while a in no or a in curr_no:
                #print('\t',a)
                if skip:
                    break;

                dm[a,ii] = 1e10
                a = np.argmin(dm[:,ii])

                if g.atom_names[a][0] == g_.atom_names[ii][0]:
                    no.append(a)
                    ope.append([a,ii])
                    u_gc_.append(gc_[ii])
                    dm[:,ii] = 1e10
                    dm[a,:] = 1e10
                    skip = True
                    break;
                else:
                    curr_no.append(a)

                cnt += 1
                #escape if something is broken
                if cnt > len(gc) + 50:
                    no.append(-1)
                    ope.append([-1,-1])
                    break;

        cgc = []
        for kk in no:
            if kk != -1:
                cgc.append(gc[kk])
        cgc = np.array(cgc)
        #print(no)
        crms = rms(cgc,np.array(u_gc_))
        if crms < r:
            r = crms
            big_no = copy.deepcopy(no)

    o = []
    for ii in big_no:
        if ii != -1:
            o.append(ii)
        #print(jj,crms)
    return r, o

def adjusted_rrms(gc,gc_):

    r = 1e+10
    nrc = []

    for jj in range(len(gc_)):
        dm = distance_matrix(gc,gc_)
        no = []
        ope = []
        curr_no = []
        for kk in range(len(gc_)):
            ii = (jj + kk) % len(gc_)
            #print(ii)

            a = np.argmin(dm[:,ii])

            skip = False
            curr_no = []

            no.append(a)
            ope.append([a,ii])
            dm[:,ii] = 1e10
            dm[a,:] = 1e10


        cgc = []
        for kk in no:
            cgc.append(gc[kk])
        cgc = np.array(cgc)
        #print(no)
        crms = rms(cgc,gc_)
        if crms < r:
            r = crms
            nrc = cgc

        #print(jj,crms)
    return r, nrc

def align_prots(i_,pr_,pc_,pr,pc):
    #get the ca of both pdbs
    y = []
    for ii in i_:
        if ii[0] not in y:
            y.append(ii[0])
    coor_ = []
    coor = []

    for ii in range(len(pr_)):
        if pr_[ii][-1] == 'CA':
            if pr_[ii][0] in y:
                coor_.append(pc_[ii])

    for ii in range(len(pr)):
        if pr[ii][-1] == 'CA':
            if pr[ii][0] in y:
                coor.append(pc[ii])

    coor = np.array(coor)
    coor_ = np.array(coor_)
    com = np.mean(coor,axis=0)
    com_ = np.mean(coor_,axis=0)

    pc_ -= com_
    pc -= com
    coor -= com
    coor_ -= com_

    H = np.dot(np.transpose(coor_),coor) / len(coor)


    U, S, V = np.linalg.svd(H)
    #print(np.shape(U),np.shape(V),np.shape(S))
    R = np.dot(U,V)

    #remove reflection
    if np.linalg.det(R) < 0:
        F = np.eye(3)
        F[2,2] = -1
        R1 = np.dot(U,F)
        R = np.dot(R1,V)

    #print(R)
    R = Rot.from_matrix(R)
    pc = R.apply(pc)
    coor = R.apply(coor)
    #print(np.mean(coor,0),np.mean(coor_,0))
    #print(coor,coor_)

    #print(rms(coor,coor_))
    return pc, pc_

def get_prot_rmsd(pc,pr,pc_,pr_):
    #get protein info
    prot_coor = []

    for ii in range(len(pr)):
        if pr[ii][-1] == 'CA':
            prot_coor.append(pc[ii])
    pc = copy.deepcopy(prot_coor)

    prot_coor = []
    for ii in range(len(pr_)):
        if pr_[ii][-1] == 'CA':
            prot_coor.append(pc_[ii])
    pc_ = copy.deepcopy(prot_coor)
    #print(len(pc),len(pc_))

    #print(pc[:5])
    #print(pc_[:5])

    return rms(pc[1:],pc_)

def fix_num(pr,pr_):
    #simplify to the simple CA only
    ca, ca_ = [], []
    for ii in pr:

        if 'CA' in ii[-1]:
            #print(ii)
            ca.append([int(ii[0]) , ii[2] ])
    for ii in pr_:

        if 'CA' in ii[-1]:
            #print('\t',ii)
            ca_.append([int(ii[0]) , ii[2] ])

    #print(ca[:10],'\n',ca_[:10])
    best_corr = 0

    best_adj = 0

    for ii in range(-500,500):
        corr = 0
        new_ca = []
        for jj in ca:
            new_ca.append([jj[0]+ii,jj[1]])

        for jj in new_ca:
            if jj in ca_:
                corr += 1

        if corr > best_corr:
            best_corr = corr
            best_adj = ii
        if corr > len(ca_) - 10:
            best_corr = corr
            best_adj = ii
            break;

    for ii in range(len(pr)):
        pr[ii][0] = str( int(pr[ii][0]) + best_adj )

    print(best_adj)

    return pr

def fnat_dice(cint,cint_):
    #print(cint)
    #print(cint_)
    c, c_ = [], []
    for ii in cint:
        for jj in ii:
            c.append(jj)
    for ii in cint_:
        for jj in ii:
            c_.append(jj)
    y = np.unique(np.array(c_).astype(int))
    y_hat = np.unique(np.array(c).astype(int))

    a = np.max(y_hat)
    if np.max(y) > a:
        a = np.max(y)

    y_arr = np.zeros(a + 500)
    y_pred_arr = np.zeros(a + 500)

    y_arr[y] = 1
    y_pred_arr[y_hat] = 1

    d = 2 * np.sum(y_arr * y_pred_arr) / (np.sum(y_arr) + np.sum(y_pred_arr))
    return d

def get_clash(pc,gc,vdw=1.75):
    dm = distance_matrix(gc,pc)
    dm = dm < vdw
    n_clash = np.sum(dm)
    return n_clash

def calc_metrics(decoy,native,same_ligand=True,is_align=True,is_same_num=True):
    """
    input:
        decoy (str): file name string of predicted structure
        native (str): file name string of native structure
        same_ligand (bool): if the ligand used is longer than the native ligand then False
    output:
        d (float): Dice of the prediction
        rirms (float): Ring RMS
        lrms (float): Ligand RMS
        dockq (float): dockq score
        s (str): string of d,rirms,lrms,dockq for easy printing
    """

    pr, pc, i, ran, rcom, gc, g, cint = get_all_info(decoy)
    pr_, pc_, i_, ran_, rcom_, gc_, g_, cint_ = get_all_info(native)

    ab_clash = get_clash(pc,gc,vdw=VDW_CUT) // 1
    aa_clash = (get_clash(gc,gc,vdw=1) - len(gc)) // 2
    #print(rcom,'\n',rcom_)
    #print(cint,'\n',cint_)

    #print(pr)

    if is_same_num == False:
        pr = fix_num(pr,pr_)

    if is_align == False:
        #print(i)
        pc, pc_ = align_prots(i_,pr_,pc_,pr,pc)
        #ca_rms = get_prot_rmsd(pc,pr,pc_,pr_)
        #print(ca_rms)



    rres = []
    rres_ = []

    #print(i,'\n\n')
    #print(i_)

    #print(len(gc),len(gc_))
    #print(len(cint),len(cint_))
    #print(len(rcom),len(rcom_))

    #print(i,i_)

    d = dice(i,i_)
    rirms = rms(rcom,rcom_)
    lrms = rms(gc,gc_)
    f = fnat(cint,cint_)

    if lrms == -1:
        same_ligand = False

    if same_ligand == False:
        #print('ppoopp')
        #change residues - look at the closest atoms to it all
        rirms, new_ring_coor = adjusted_rrms(rcom,rcom_)
        #cint = find_interactRingRes(new_ring_coor,pc,pr,INTERACT = INTERACT_RING)
        #print('r',cint)
        cint = find_interactRingAtomRes(new_ring_coor,g,pc,pr,INTERACT = INTERACT_RING)
        #print('a',cint,'\n',cint_)

        f = adjusted_fnat(cint,cint_)

        lrms, atoms = adjusted_lrms(gc,gc_,g,g_)
        #print(len(gc),len(gc_),len(atoms))
        #print('\t\t',lrms)
        #print(rirms,lrms)
        new_gc = []
        for ii in range(len(atoms)):
            new_gc.append(gc[atoms[ii]])
        new_gc = np.array(new_gc)
        #print(atoms)

        #recalculate the dice based on info
        i = find_interactChains(new_gc,pc,pr,pr)
        #print(i)
        d = fnat_dice(cint,cint_)



    #s = str(round(d,3)) + ',' + str(round(rirms,3)) + ',' + str(round(lrms,3)) + ',' + str(round(dockq,3)) + '\n'

    return d,f,rirms,lrms,ab_clash,aa_clash
