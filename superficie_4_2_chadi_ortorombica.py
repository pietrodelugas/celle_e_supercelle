# %%
import pymatgen.core as pmg
import numpy as np

# %% [markdown]
# First thing we generate the bulk $\mathrm{fcc}$ structure. 
# * we define the lattice with the 3 fcc vectors
# * a list with a symbol for the atoms in the structure ( in this just one Al) 
# * the coordinates in the same order as the symbols is a list of np.arrays of length 3

# %%
alat = 10.21 
lattice = np.array([[0.,1.,1.], [1.,0.,1], [1.,1.,0]])* alat *0.5
species = ['Si','Si']
coords=[np.array([0,0,0]), np.array([0.25,0.25,0.25])]

s = pmg.Structure(lattice=lattice, species=species , coords=coords, coords_are_cartesian=False) 

# %% [markdown]
# just to check everything is fine we make print out the info about the structure

# %%
s

# %% [markdown]
# we now generate a supercell with the $z$  axis oriented along the $[1 0 0]$ direction and the other 2 axes perpendicular to it, generating the 4X2 conventional cell  
# * we copy the bulk structure into s1 
# * we transform `sup1` in a supercell. The  input in the second line means that the new vectors of the lattice for the supercell are $$\bf{a_1^\prime} = -\bf{a_1} + -2\bf{a_2}+2\bf{a_3}  \; ; \; \bf{a_2^\prime} = \bf{a_1}-2\bf{a_2}+2\bf{a_3} \: ; \: \bf{a_3^\prime} = -\bf{a_1} + \bf{a_2} +\bf{a_3}  $$

# %%
sup1 = s.copy()
sup1.make_supercell([[-1,-2,2],[1,-2,2],[-1,1,1]], to_unit_cell = True)

# %% [markdown]
# Now we need to transform the axes in order that $\bf{c}$ points in the $z$ direction.
# The object `lattice` of pymatgen has all the info that we need. The 3 axis lengths $a,b,c$ and the angles between the axes $\alpha$ between b and c axes , $\beta$ between $\bf{a}$  and $\bf{c}$ and $\gamma$ between $\bf{a}$ and $\bf{b}$. 
# 
# We just set ${\bf c}^\prime = c \cdot \hat{\bf z}$; we then put ${\bf a}^\prime = a\cdot sin(\beta)\cdot \hat{\bf x} + a\cdot cos(\beta)\cdot\hat{\bf z}$ in the $xz$ plane forming an angle $\beta$ with ${\bf c}^\prime$; the components of $\mathbf{b}^\prime$ are found imposing that ${\bf b^\prime \cdot c^\prime} = b\,c\cdot cos(\alpha)$ and 
# ${\bf b^\prime \cdot a^\prime} = a\,b \cdot cos(\gamma)$
# 
# Just use this numbers to create the new vectors and print them to check we have done right

# %%
latt = sup1.lattice
a = latt.a ; alpha = latt.alpha*np.pi/180.0
b = latt.b ; beta = latt.beta*np.pi/180.0
c = latt.c ; gamma = latt.gamma *np.pi/180.0
print(a,b,c, alpha/np.pi*180, beta/np.pi*180, gamma/np.pi*180)
newa3 = c * np.array([0,0,1])
newa1 = a * np.array([np.sin(beta),0,np.cos(beta)])
ca = (np.cos(gamma) - np.cos(beta)*np.cos(alpha))/np.sin(beta)
cb = np.sqrt(np.sin(alpha)**2-ca**2)
newa2 = b * np.array([ca,cb,np.cos(alpha)])

# %% [markdown]
# 
# we rotate the new $\bf{a_1}$ and $\bf{a_2}$ in the plase so that these  axes are oriented compatibly with `pw.x` for the orthorhombic base-centered case with `ibrav=-9` 
# 
# ``` 
# (a/2,-b/2,0),(a/2,b/2,0),(0,0,c)
# ```

# %%
bisect = (newa1 + newa2)/np.sqrt((newa1+newa2).dot(newa1+newa2)) 
cbisect = newa1.dot(bisect) / a 
newa1 = a * np.array([cbisect, - np.sqrt(1 - cbisect * cbisect), 0.0]) 
newa2 = a * np.array([cbisect,   np.sqrt(1 - cbisect * cbisect), 0.0])
newA = 4 * np.sqrt(1./2) * alat 
newB = 2 * np.sqrt(1./2) * alat 

# %% [markdown]
# now we replace these three axes in the structure.  then we print out structure info  to compare lengths, angles and fractional coordinates. 

# %%
sup1.lattice = pmg.Lattice(np.array([newa1, newa2,newa3]))
sup1, sup1.lattice

# %% [markdown]
# Now that the c vector is  oriented along the z cartesian axis we can build any supercell along the $[100]$ direction of any the periodicity along the xy  plane and of any thickess along the z direction 
# 
# the first two vectors define the periodicity in the plane , the third vector the thickness. For example a $2\times2$ surface with thickness 3 is generated like this:
# 
# * first we generate the supercell 

# %%
sup2 = sup1.copy()
sup2.make_supercell([[1,0,0],[0,1,0],[0,0,8]], to_unit_cell = False)
sup2

# %% [markdown]
# * then we extract the coordinates and order them along the z direction, in case we have different species we need to sort also the symbols .... 

# %%
coords = sup2.cart_coords
species = sup2.species
temp  = sorted(zip(species,coords), key= lambda x: x[1][2])
species = [_[0] for _ in temp ]
coords  = [_[1] for _ in temp ]
coords

# %% [markdown]
# * then we modify the periodicit along the a3 axes. We replace previous $\bf{c}$ vector with one that is longer than the lattice periodity in that direction this create a gap between the slabs. 

# %%
a1 = sup2.lattice.matrix[0]
a2 = sup2.lattice.matrix[1]
vacuum = 20 # length of the vacuum region 
a3 = sup2.lattice.matrix[2]+np.array([0.,0., vacuum])
surface = pmg.Structure(lattice = np.array([a1,a2,a3]), coords =coords, species= species, coords_are_cartesian=True)
from	collections import deque
it = (f"{_.specie.name}   {_.coords[0]:12.6f}  {_.coords[1]:12.6f}  {_.coords[2]:12.6f} \n" for _ in surface.sites) 
celldm1 = newA 
celldm2 = newB/newA 
celldm3 = np.sqrt(surface.lattice.matrix[2].dot(surface.lattice.matrix[2]))/celldm1  
with open('posizioniSi100_4X2_chadi_orthorombica','w') as f:
	f.write(f"ibrav=-9    celldm(1)={celldm1:9.6f}   celldm(2)={celldm2:9.6f}   celldm(3)={celldm3:9.6f}" +
	        f"   nat={len(surface.sites)}\n\n")
	s = [" ".join((f"{round(_,6):12.6}" for _ in surface.lattice.matrix[i])) for i in [0,1,2]]    
	f.write ( "\n".join(s) + "\n\n")
	deque((f.write(_) for _ in it )) 
 

# %%
surface.to(filename='surf.cif')

# %%
p = surface.sites[4]

# %%
surface.lattice.matrix

# %%
a = surface.lattice.matrix 
print (f"alat = { np.sqrt(a[0].dot(a[0]))} cosgamma={a[0].dot(a[1])/a[1].dot(a[1])} ")

# %%
np.arccos(0.6)/np.pi**2

# %%



