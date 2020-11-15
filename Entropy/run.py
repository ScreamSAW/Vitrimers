import math
import hoomd
from hoomd import md
from hoomd import group
from hoomd import deprecated
from hoomd import data
from hoomd import analyze
from hoomd import dump

#initialize
hoomd.context.initialize();

#box length and numbers
length_box=float(40)
NPoly=100
NXlinker=50

#compute the box volume
varSysVol = length_box * length_box * length_box

#steps
totrun=1000

#define the polymer bonds
PolyBondStruct = []
	#define the backbone
PolyBondStruct.append((0, 1, 'real'))      # 0 is the starter node
PolyBondStruct.append((1, 2, 'real'))
PolyBondStruct.append((2, 3, 'real'))
PolyBondStruct.append((3, 4, 'real'))
PolyBondStruct.append((4, 5, 'real'))
PolyBondStruct.append((5, 6, 'real'))
PolyBondStruct.append((6, 7, 'real'))
PolyBondStruct.append((7, 8, 'real'))
	#define the x-link points
PolyBondStruct.append((1, 9, 'real'))
PolyBondStruct.append((4, 10, 'real'))
PolyBondStruct.append((7, 11, 'real'))
	#defien the -BA structure at x-link points
PolyBondStruct.append((9, 12, 'fake'))
PolyBondStruct.append((10, 13, 'fake'))
PolyBondStruct.append((11, 14, 'fake'))

#define the x-linker
XBondStruct = []
XBondStruct.append((0,1, 'real'))
XBondStruct.append((1,2, 'real'))

#build molecules
polymerA = dict(bond_len=2.1, type=(9*['N']+3*['B']+3*['A']), bond=PolyBondStruct, count=NPoly)
Xlinker = dict(bond_len=2.1, type=(['C']+['N']+['C']), bond=XBondStruct, count=NXlinker)

#generate the polymer system
system = deprecated.init.create_random_polymers(box=data.boxdim(volume=varSysVol), polymers=[polymerA,Xlinker], separation=dict(B=1, N=1, A=1, C=1),seed=9527)

#neighboorlist
nl=md.nlist.cell()
all=group.all()

#RevCross 3-body potential
rc=md.pair.revcross(r_cut=1.5,nlist=nl)
rc.pair_coeff.set(['N','A','B','C'],['N','A','B','C'],sigma=0,n=0,epsilon=0,lambda3=0)
rc.pair_coeff.set('A','B',sigma=0.5,n=10,epsilon=100,lambda3=1)
rc.pair_coeff.set('C','B',sigma=0.5,n=10,epsilon=100,lambda3=1)

#rigid


#L-J 2-body potential
sigmalj=0.9
my_rcut=math.pow(2,1/6)*sigmalj
potRep=md.pair.lj(r_cut=my_rcut,nlist=nl)
potRep.pair_coeff.set(['N','A','B','C'],['N','A','B','C'],epsilon=1.,sigma=sigmalj)
potRep.pair_coeff.set('A','B',epsilon=0.0,sigma=0.01)
potRep.pair_coeff.set('C','B',epsilon=0.0,sigma=0.01)
potRep.set_params(mode="shift")

#bond coeff
harm=md.bond.harmonic()
harm.bond_coeff.set('real',k=1000,r0=1)
harm.bond_coeff.set('fake',k=1000,r0=1)

#output file
NAME='50.log'
og1 = analyze.log(filename=NAME, quantities=['time', 'pressure', 'potential_energy', 'temperature'], period=1000, overwrite=True)

#temperature
T=298.0


#dump file
deprecated.dump.xml(group=all, filename="dump.xml", position=True, velocity=True, type=True, period=None, restart=False)
dump.dcd(filename="movie.dcd",period=1000)

'''
#### DUMP FOR VMD
deprecated.dump.xml(group=all, filename="dump.xml", position=True, velocity=True, type=True, period=None, restart=False)
dump.dcd(filename="movie.dcd",period=1000)
'''
##FIRST RUN WITH LIMIT TO EQUILIBRATE
md.integrate.mode_standard(dt=0.0001)
a=md.integrate.nve(group=all, limit=0.001)
hoomd.run(int(totrun)*10)
md.integrate.nve.disable(a)
#free A in -BA
harm.bond_coeff.set('fake',k=0,r0=0)

# output the snapshot
deprecated.dump.xml(group=group.all(), filename="_Init.xml", vis=True)
NAME='50.gsd'
dump.gsd(filename=NAME, group=group.all(), truncate=True, period=100, phase=0,overwrite=True)

#temperature quench NVT
md.integrate.mode_standard(dt=0.0001)
b=md.integrate.nvt(group=all, kT=T, tau=1) 
hoomd.run(int(totrun)*100000)
md.integrate.nvt.disable(b)

