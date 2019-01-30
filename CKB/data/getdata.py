class data():
	"""data?"""
	def __init__(self, serial, angle, x, y):
		# super(data, self).__init__()
		self.serial = serial
		self.angle = angle
		self.x = x
		self.y = y
		self.ldr = y/x
		self.force = (x**2+y**2)**.5
		
f = open('force.dat','r')
data_list = f.read().split('\n')
# for i in range(len(data_list)):
i=0
datas = []
while True:
	if data_list[i] == '':
		break
	a = data_list[i].split(' ')
	i+=1
	b = data_list[i].split(' ')
	i+=1
	datas.append(data(serial=int(a[0]),angle=int(a[1]),x=float(b[0]),y=float(b[1])))

tmp = 0
for d in datas:
	print d.ldr
	if tmp<d.ldr:
		tmp = d.ldr
		serial = d.serial
		ang = d.angle
print 'biggest lift drag ratio = ', tmp, serial, ang
# print data_list