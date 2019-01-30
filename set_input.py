import logging
logging.basicConfig(level=logging.DEBUG,
                    filename ='autosimulate.log',
                    format='%(asctime)s - %(levelname)s : %(message)s')
# Objective:
# get AOA vs lift & drag data for every airfoil
# save ib & pressure for visualization
#
# Algorithm:
#
# Generate NACA .inp files: CKB/NACA/foil/NACAxxxx.inp
# for each NACAxxxx:
# 	make cleandata
# 	copy NACAxxxx.inp to input/body.001.inp
# 	set AOA in input/ib.inp to -5
# 	while AOA<=25:
# 		execute bin/ib
# 		move output/snapshots/ib0001000.var to CKB/data/ib/[NACA]_[AOA].var
# 		move output/snapshots/pressure0001000.var to CKB/data/pressure/[NACA]_[AOA].var
# 		read output/responds/simple_force.dat, record last line data
# 		AOA += 5
# 	end while
# 	create file: CKB/data contains lift and drag for each AOA
# end for
import os
import sys



# predefined constant
current = os.getcwd()
# print(current) # for debuging
foil = '/CKB/NACA/foil/'
config = """ &READ_PARAMETERS
 ISTART = 0,
 ISTOP = 1000,
 ISAVE = 100,
 DT = 0.005,
 M = 200,
 N = 200,
 LEN = 3,
 OFFSETX = 1,
 OFFSETY = 1.5,
 MGRIDLEV = 5,
 RE = 100,
 NCR = 1,
 OUTPUT_SF = T,
 COMPUTE_PRESSURE = T,
 AOA = %d
/
"""
# get generated NACA list
naca_log = current + foil + 'naca_gen.log'
naca_log_file = open(naca_log, 'r')

num_list = []
for line in naca_log_file:
	num_list.append(line[5:9])
naca_log_file.close()
# print num_list
#num_list = [num_list[0]] # for debuging

# in log file: 1st line: NACA serial
#			   2nd line: AOA
main_log_file = open(current+'/main_log.log','r')
s_mainlog = main_log_file.read()
lis_mainlog = s_mainlog.split('\n')

num = lis_mainlog[0]
AOA = int(lis_mainlog[1])
# print 'Writing body.001.inp for NACA %s, AOA = %d'%(num,AOA)

# os.system('make cleandata')
foilinp = current + foil + 'NACA' + num + '.inp'
command = 'cp '+ foilinp + ' ' + current + '/input/body.001.inp'
logging.info('copying %s'% num)
os.system(command)
#AOA = 25

#print(config%AOA)
# print 'Writing ib.inp for NACA %s, AOA = %d'%(num,AOA)
config_file = open(current+'/input/ib.inp', 'w')
config_file.write(config % AOA)
config_file.close()