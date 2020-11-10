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



length_box=float(40)
NAstar=950
NBstar=550


# compute the box volume
varSysVol = length_box * length_box * length_box
totrun=10000

#define the star structure
varBondStruct = []
varBondStruct.append((0, 1, 'polymerCha'))      # 0 is the central node
varBondStruct.append((0, 4, 'polymerCha'))
varBondStruct.append((0, 7, 'polymerCha'))
varBondStruct.append((0, 10, 'polymerCha'))
varBondStruct.append((0, 13, 'polymerCha'))
varBondStruct.append((0, 16, 'polymerCha'))
varBondStruct.append((0, 19, 'polymerCha'))
varBondStruct.append((0, 22, 'polymerCha'))
#8 arms
varBondStruct.append((1, 2, 'polymerCha'))
varBondStruct.append((2, 3, 'polymerCha'))

varBondStruct.append((4, 5, 'polymerCha'))
varBondStruct.append((5, 6, 'polymerCha'))

varBondStruct.append((7, 8, 'polymerCha'))
varBondStruct.append((8, 9, 'polymerCha'))

varBondStruct.append((11, 12, 'polymerCha'))
varBondStruct.append((10, 11, 'polymerCha'))

varBondStruct.append((13, 14, 'polymerCha'))
varBondStruct.append((14, 15, 'polymerCha'))

varBondStruct.append((16, 17, 'polymerCha'))
varBondStruct.append((17, 18, 'polymerCha'))

varBondStruct.append((19, 20, 'polymerCha'))
varBondStruct.append((20, 21, 'polymerCha'))

varBondStruct.append((22, 23, 'polymerCha'))
varBondStruct.append((23, 24, 'polymerCha'))

polymerStarA = dict(bond_len=0.5, type=(['C'] + 2*['N'] + ['A'] + 2*['N'] + ['A'] + 2*['N'] + ['A']+ 2*['N'] + ['A'] + 2*['N'] + ['A'] + 2*['N'] + ['A'] + 2*['N'] + ['A'] + 2*['N'] + ['B']), bond=varBondStruct, count=NAstar)
polymerStarB = dict(bond_len=0.5, type=(['C'] + 2*['N'] + ['B'] + 2*['N'] + ['B'] + 2*['N'] + ['B']+ 2*['N'] + ['B'] + 2*['N'] + ['B'] + 2*['N'] + ['B'] + 2*['N'] + ['B'] + 2*['N'] + ['A']),  bond=varBondStruct, count=NBstar)


#generate the polymer system
system = deprecated.init.create_random_polymers(box=data.boxdim(volume=varSysVol), polymers=[polymerStarA,polymerStarB], separation=dict(C=0.2,B=0.2, N=0.2, A=0.2),seed=9527)

#neighboorlist
nl=md.nlist.cell()
all=group.all()

#RevCross 3-body potential
rc=md.pair.revcross(r_cut=1.8,nlist=nl)
rc.pair_coeff.set(['N','A','B','C'],['N','A','B','C'],sigma=0,n=0,epsilon=0,lambda3=0)
rc.pair_coeff.set('A','B',sigma=0.5,n=10,epsilon=100,lambda3=1)

#L-J 2-body potential
sigmalj=0.9
my_rcut=math.pow(2,1/6)*sigmalj
potRep=md.pair.lj(r_cut=my_rcut,nlist=nl)
potRep.pair_coeff.set(['N','A','B','C'],['N','A','B','C'],epsilon=1.,sigma=sigmalj)
potRep.pair_coeff.set('A','B',epsilon=0.0,sigma=0.01)
potRep.set_params(mode="shift")

harm=md.bond.harmonic()
harm.bond_coeff.set('polymerCha',k=1000,r0=1)

#output file
NAME='10-2.log'
og1 = analyze.log(filename=NAME, quantities=['time', 'pressure_xy','pressure_xz','pressure_yz'], period=1, overwrite=True)

#temperature
T=298.0

'''
#dump file
deprecated.dump.xml(group=all, filename="dump.xml", position=True, velocity=True, type=True, period=None, restart=False)
dump.dcd(filename="movie.dcd",period=1000)
'''
'''
#### DUMP FOR VMD
deprecated.dump.xml(group=all, filename="dump.xml", position=True, velocity=True, type=True, period=None, restart=False)
dump.dcd(filename="movie.dcd",period=1000)
'''
##FIRST RUN WITH LIMIT TO EQUILIBRATE
md.integrate.mode_standard(dt=0.0001)
a=md.integrate.nve(group=all, limit=0.001)
hoomd.run(int(totrun))

md.integrate.nve.disable(a)

# output the snapshot
deprecated.dump.xml(group=group.all(), filename="_Init.xml", vis=True)
#NAME='%scl/restart%s.gsd'%(sys.argv[2],sys.argv[1])
NAME='10.gsd'
dump.gsd(filename=NAME, group=group.all(), truncate=True, period=100000, phase=0,overwrite=True)


##temperature quench NVT
md.integrate.mode_standard(dt=0.0001)
b=md.integrate.nvt(group=all, kT=T, tau=1) 
hoomd.run(int(totrun)*5000)
md.integrate.nvt.disable(b)

