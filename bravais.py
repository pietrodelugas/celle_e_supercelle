"""Bravais module

contains many small functions use to manipulate Bravais lattices,
and atomi structure info in general
"""
import numpy as np
# funzioncine ricorrenti: 
vol   = lambda a: a[0].dot(np.cross(a[1],a[2]))
#
norm  = lambda v: np.sqrt(v.dot(v))
#
recip = lambda a: 2.0 * np.pi /vol(a) * np.array([np.cross(a[i-1],a[j-1]) for i,j in [(2,3),(3,1),(1,2)]])
#
nm    = lambda at_sc,at_0,d: tuple(
            (s*(int(round(norm(recip(at_0)[d].dot(at_sc)/2/np.pi),0))+2) for s in (-1,1))
)
#
t0_sc_frac = lambda at_0, at_sc, i,j,k: recip(at_sc).dot( np.array([i,j,k]).dot(at_0))/2./np.pi
#
inbox =      lambda at_0, at_sc,i,j,k: (
    (t0_sc_frac(at_0,at_sc,i,j,k) > -1.e-6   ).all() and 
    (t0_sc_frac(at_0,at_sc,i,j,k) < 1.0-1.e-6).all()
) 
#
trasl = lambda at_sc,at_0: filter(lambda t: inbox(at_0,at_sc,*t), ((i,j,k) 
                                        for i in range(*nm(at_sc,at_0,0)) 
                                        for j in range (*nm(at_sc,at_0,1))
                                        for k in range(*nm(at_sc, at_0,2))
                                       ))
#
frac_pos = lambda at, pos:  [(_[0],recip(at).dot(_[1])/2./np.pi) for _ in pos]
#
scale_pos = lambda fac, pos: [(_[0],_[1]/fac) for _ in pos]
#
posline = lambda v: " ".join([f"{x_:12.7f}" for x_ in v])
write_pos_iteration =  lambda pos,f: (f.write(f"{_[0]:4s}" + posline(_[1])+"\n") for _ in pos)
#
tetragonal_lattice = lambda alat,bsua, csua: alat*np.array([[1.,0.,0.],[0.,bsua,0.],[0.,0.,csua]])
hexagonal_lattice  = lambda alat,csua: alat*np.array([[1.,0.,0.],[-0.5, np.sqrt(0.75),0.],[0., 0.,csua]]) 
#
rotate_pos  = lambda old, new, pos: [(_[0],_[1].dot(new)) for _ in frac_pos(old,pos)] 
#
def supercell_pos(at_sc, at_0, pos0):
    pos_frac0 = frac_pos(at_0, pos0) 
    return  [(_[0],(_[1]+np.array(t)).dot(at_0)) for _ in pos_frac0 for t in trasl(at_sc, at_0) ]
#
def cubic_ibrav_alat(at):
    from sys import stdout
    #
    if round(at[0].dot(at[1]),8) == 0:
        alat = norm(at[0])
        ibrav = 1
    elif round(at[0].dot(at[1])/norm(at[0])/norm(at[1]),8) == 0.5:
        alat = norm(at[0])*np.sqrt(2.)
        ibrav = 2
    elif round(at[0].dot(at[1])/norm(at[0])/norm(at[1]),8) == round(-1./3.,8):
        alat = norm(at[0])/0.8660254037844386
        ibrav = -3 
    return alat, ibrav 


def tetragonal_ibrav_alat(at):
    alat = norm(at[0])
    bsua = norm(at[1])/alat
    csua = norm(at[2])/alat 
    if round(bsua,8) == 1.0:
        ibrav = 6 
    else:
        ibrav = 8 
    return alat, bsua, csua, ibrav 
#
def hexagonal_ibrav_alat(at): 
  alat = norm(at[0])
  csua = norm(at[2])/norm(at[0])
  try:
    assert ( round(at[0].dot(at[1])/norm(at[0])**2,8) == -0.5 and  
             round(norm(np.cross(at[0],at[1]))/norm(at[0])**2,8) == round(np.sqrt(0.75),8)
           )
  except AssertionError:
    print (at[0].dot(at[1])/norm(at[0])**2 , norm(np.cross(at[0],at[1]))/norm(at[0])**2)
    raise AssertionError
  return  alat, csua, 6 
    
#
def cubic_alat_pos(at, pos,filename=None):
    from collections import deque
    from sys import stdout 
    #
    alat, ibrav = cubic_ibrav_alat(at)
    #
    scaled_pos = scale_pos(alat, pos)
    if filename is None:
        print (f"alat = {alat:10.6f}")
        print (f"ibrav = {ibrav:4d}" )
        print (f"nat = {len(pos):4d}\n") 
        deque(write_pos_iteration(scaled_pos, stdout), maxlen=0)
        return 1
    else:
        with open(filename,'w') as f:
            f.write(f"alat = {alat:10.6f}\n")
            f.write(f"ibrav = {ibrav:4d}\n")
            f.write(f"nat = {len(pos):4d}\n\n")
            deque(write_pos_iteration(scaled_pos,f),maxlen=0)
        return 1 
#    
def tetragonal_alat_pos(at, pos, filename=None):
    from collections import deque
    from sys import stdout
    #
    alat, bsua, csua, ibrav = tetragonal_ibrav_alat(at)
    #
    scaled_pos = scale_pos(alat,pos)
    if filename is None:
        print (f"alat = {alat:10.6f}")
        if ibrav == 8: print(f"bsua = {bsua:9.6f}")
        print (f"csua = {csua:9.6f}")
        print (f"ibrav = {ibrav:4d}")
        print (f"nat = {len(pos):4d}\n")
        deque(write_pos_iteration(scaled_pos, stdout))
        return 1
    else:
        with open(filename,'w') as f:
            f.write (f"alat = {alat:10.6f}\n")
            if ibrav == 8: f.write(f"bsua = {bsua:9.6f}\n")
            f.write(f"csua = {csua:9.6f}\n")
            f.write(f"ibrav = {ibrav:4d}\n")
            f.write(f"nat = {len(pos):4d}\n\n")
            deque(write_pos_iteration(scaled_pos, f))
            return 1
    
    
                
        
     