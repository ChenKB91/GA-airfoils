import os
import sys
import logging
logging.basicConfig(level=logging.DEBUG,
                    filename ='autosimulate.log',
                    format='%(asctime)s - %(levelname)s : %(message)s')
current = os.getcwd()
sys.__stderzr__ = None

serials = ['0012', '0020', '0030', '0050', 
		   '2212', '2220', '2230', '2250', 
		   '2412', '2420', '2430', '2450', 
		   '2612', '2620', '2630', '2650', 
		   '4212', '4220', '4230', '4250', 
		   '4412', '4420', '4430', '4450', 
		   '4612', '4620', '4630', '4650', 
		   '6212', '6220', '6230', '6250', 
		   '6412', '6420', '6430', '6450', 
		   '6612', '6620', '6630', '6650']

data = open(current+'/output/responds/simple_force.dat','r')
lis = data.read().split('\n')
data.close()
s = lis[999]
l1 = [x for x in s.split(' ') if x != '']

# in log file: 1st line: NACA serial
#			   2nd line: AOA
main_log_file = open(current+'/main_log.log','r')
# s_mainlog = main_log_file.read()
lis_mainlog = main_log_file.read().split('\n')
main_log_file.close()

num = lis_mainlog[0]
AOA = int(lis_mainlog[1])

logging.info('collect output: writing force.dat for %s %d\n'%(num,AOA))

os.system('cp %s/output/snapshots/ib0000900.var %s/CKB/data/ib/%s_%s.var' % (current,current,num,AOA))
os.system('cp %s/output/snapshots/pressure0000900.var %s/CKB/data/pressure/%s_%s.var' % (current,current,num, AOA))


# print 'force: %s %s'%(l1[1],l1[2])
forcefile = open(current+'/CKB/data/force.dat','a')
forcefile.write(num+' '+str(AOA)+'\n')
forcefile.write(l1[1]+' '+l1[2]+'\n')
forcefile.close()

brk = False
if AOA < 25:
	AOA+=5
elif AOA == 25:
	i = serials.index(num)+1
	AOA = -5
	if i >= len(serials):
		logging.warning('STOP!!!')
		brk = True
		# sys.__stderr__ = "999"
	else:
		num = serials[i]
else:
	logging.error('unknown error')

main_log_file2 = open(current+'/main_log.log','w')
main_log_file2.write(num+'\n')
main_log_file2.write(str(AOA))
main_log_file2.close()

sys.stderr.write('0')
if brk:
	# print 999
	sys.stderr.write('999')


