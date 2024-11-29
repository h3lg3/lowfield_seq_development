import pypulseq as pp

from packages.mr_systems import lumina as system

seq = pp.Sequence(system=system)
seq.read(r'E:\Python\tse_3d_lumina_64_64_64.seq', detect_rf_use = True)




print(seq.definitions)