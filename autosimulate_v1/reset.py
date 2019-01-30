import os

current = os.getcwd()
mlog = open(current+'/main_log.log','w')
mlog.write('0012\n-5')
mlog.close()
open(current+'/CKB/data/force.dat','w').close()
open(current+'/autosimulate.log','w').close()
f=open(current+'/input/ib.inp','w')
s=''' &READ_PARAMETERS
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
 AOA = -5
/
'''
f.write(s)
f.close()