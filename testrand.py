import random as rd

n = 120000
t = 0
for i in range(n):
	a = rd.randint(1,6)
	b = rd.randint(1,6)
	c = rd.randint(1,6)
	if a!=b and b!=c and a!=c:
		t+=1

print t