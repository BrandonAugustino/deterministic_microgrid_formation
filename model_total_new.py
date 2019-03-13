from sympy import binomial
from gurobipy import *
m = Model("thetamodelfull")

#######    DATA LOAD    ###############
# file = open("/Users/mertcanyetkin/Dropbox/SmartCities/microgrids-code/resil/ieee30_s5_python.dat","r")
file = open("../resil/ieee30_s5_python.dat","r")


data = file.read().split('\n')
N = int(data[0].split()[2])                        # number of nodes
K = int(data[1].split()[2])                        # number of generators
K = [str(i+1) for i in range(K)]
BB = list(map(str, data[2].split()[2:]))           # list of nodes in network
B = set(BB)                                        # set of nodes in network
L = binomial(N,2)                                  # number of possible links

print ('test1')
Alist = []
A = []
for i in range(3,4):
	AAA1 = data[i].split()[2]                      # set of all the links in scenario i
	AA1 = AAA1.split(';')
	A1 = set(AA1)
	Alist1 = {}
	for k in BB:
		Atemp = []	
		for j in AA1:
			if list(map(int,j.split(',')))[0] == int(k):
				Atemp.append(list(map(str,j.split(',')))[1])
		Alist1[k] = Atemp
	Alist.append(Alist1)
	A.append(AA1)


AB =[]
for i in range(len(A[0])):
	AB.append(tuple(A[0][i].split(',')))


print ('test2')
p = list(map(float, data[4].split()[2:]))         # active load of each node
q = list(map(float, data[5].split()[2:]))         # reactive load of each node
w = list(map(float, data[6].split()[2:]))         # reactive load of each node --- Should this be criticality weight?
rtemp = list(map(float, data[7].split()[2:]))     # line resistance on each link
r = {}
for p1 in range(0,int(len(rtemp)/3)):
	a = (str(int(rtemp[3*p1])), str(int(rtemp[3*p1+1])))
	r[a] = rtemp[3*p1+2]


xtemp = list(map(float,data[8].split()[2:]))      # line reactance on each link -- reactance or resistance?
x = {}
for p1 in range(0,int(len(xtemp)/3)):
	a = (str(int(xtemp[3*p1])), str(int(xtemp[3*p1+1])))
	x[a] = xtemp[3*p1+2]


print ('test4')
Vo = int(data[9].split()[2])                      # Vo
VR = int(data[10].split()[2])                     # VR
epslison = float(data[11].split()[2])             # epslison
Pmax = list(map(float, data[12].split()[2:]))     # maximum generation of each generator
Qmax = list(map(float, data[13].split()[2:]))     # maximum generation of each generator
TP = int(data[14].split()[2])                     # TP
TQ = int(data[15].split()[2])                     # TQ


#dd = list(map(float, data[17].split()[2:]))      # reactive load of each node
dd = [1 for i in range(N)]  
hp = 0.1
nodeo = list(map(str, data[18].split()[2:]))
nodec = list(map(str, data[19].split()[2:]))

linko = []
AAA1 = data[20].split()[2]                       # set of all the links in scenario i
AA1 = AAA1.split(';')
linko.append(AA1)

linkoAB =[]
for i in range(len(linko[0])):
	linkoAB.append(tuple(linko[0][i].split(',')))


linkc = []
AAA1 = data[21].split()[2]                      # set of all the links in scenario i
AA1 = AAA1.split(';')
linkc.append(AA1)

linkcAB =[]
for i in range(len(linkc[0])):
	linkcAB.append(tuple(linkc[0][i].split(',')))


genlocdata = data[22].split()[2].split(';')
genloc = {}
for item in genlocdata:
	itemtemp = item.split(',')
	genloc[itemtemp[0]] = itemtemp[1]


print ('test5')
##########   Decision Variables   #############

s = m.addVars(BB, vtype=GRB.BINARY,name = "load")
z = m.addVars(K,BB,vtype = GRB.BINARY, name = "location")
b = m.addVars(AB,vtype = GRB.BINARY,name = "link")
Pg = m.addVars(BB, vtype = GRB.CONTINUOUS, name = "activepowergen")
P = m.addVars(AB, lb = -GRB.INFINITY, ub = GRB.INFINITY, vtype = GRB.CONTINUOUS, name = "activepower")
Qg = m.addVars(BB,vtype = GRB.CONTINUOUS, name = "reactivepowergen")
Q = m.addVars(AB,lb = -GRB.INFINITY, ub = GRB.INFINITY,vtype = GRB.CONTINUOUS, name = "reactivepower")
V = m.addVars(BB, lb = -GRB.INFINITY, vtype = GRB.CONTINUOUS, name = "voltage")
delta = m.addVars(AB, lb = -GRB.INFINITY, vtype = GRB.CONTINUOUS, name = "slack")
deltal = m.addVars(AB, lb = -GRB.INFINITY, vtype = GRB.CONTINUOUS, name = "slack1")

###########################################################################################
###########################################################################################
print ('test6')



############ for microgrid(2)
v = m.addVars(BB,K, vtype=GRB.BINARY,name = "node_to_microgrid")
c = m.addVars(AB,K, vtype=GRB.BINARY,name = "link_to_microgrid")

# ###########################################################################################
# ###########################################################################################
# ############################### For dispatchable

pd = m.addVars(BB,vtype = GRB.CONTINUOUS, name = "dispatch")
ps = m.addVars(BB,vtype = GRB.CONTINUOUS, name = "servedload")
qs = m.addVars(BB,vtype = GRB.CONTINUOUS, name = "servedloadq")

###########################################################################




############  Constraints     #############
######### node P constraint
m.addConstrs(((quicksum(P[i,j] for j in Alist[0][i])) ==  -ps[i] + Pg[i] for i in BB  ),"nodeP")  
######### node Q constraint
m.addConstrs(((quicksum(Q[i,j] for j in Alist[0][i])) ==  -qs[i] +  Qg[i] for i in BB ),"nodeQ")  

######### LINE CONSTRAINT  ############# 
for i,j in AB:
	m.addConstr(-TP*b[i,j] <= P[i,j])
	m.addConstr(P[i,j] <= TP*b[i,j] )
	m.addConstr(-TQ*b[i,j] <= Q[i,j] )
	m.addConstr(Q[i,j] <= TQ*b[i,j] )
	m.addConstr(P[i,j] == - P[j,i])
	m.addConstr(Q[i,j] == - Q[j,i])


print ('test7')

######### confinement of generators
m.addConstrs(0 <= Pg[j] for j in BB )
m.addConstrs(Pg[j] <= quicksum((z[i,j]*Pmax[int(i)-1]) for i in K) for j in BB )
m.addConstrs(0 <= Qg[j]  for j in BB )
m.addConstrs(Qg[j] <= quicksum((z[i,j]*Qmax[int(i)-1]) for i in K)  for j in BB )
#print(BB)

# ########### CONSTRAINTS OF VOLTAGE
for i,j in AB:
	m.addConstr(V[i] == V[j] + (r[(i,j)]*P[i,j] + x[(i,j)]*Q[i,j])/Vo + delta[i,j] + deltal[i,j]) 
	m.addConstr((-1 + b[i,j])*Vo <= delta[i,j])
	m.addConstr( delta[i,j] <= (1 - b[i,j])*Vo)
	m.addConstr(deltal[i,j] >= -0.01)
	m.addConstr(deltal[i,j] <= 0.01)

# # # ###############################################3

m.addConstrs(Vo*quicksum(z[i,j] for i in K) <= V[j]  for j in BB )
m.addConstrs( V[j] <= Vo for j in BB )
print ('test8')

# #  ###############################################3
m.addConstrs( (1 - epslison)*VR <= V[j] for j in BB )
m.addConstrs(  V[j] <= (1 + epslison)*VR  for j in BB )

############ LOGICAL CONSTRAINTS
m.addConstrs(quicksum(z[i,j] for j in BB) == 1 for i in K)

##############new 
m.addConstrs(quicksum(z[i,j] for i in K) <= 1 for j in BB)

# m.addConstr(Pg[0] == 40)
print (Pmax)



for i in genloc.keys():
	m.addConstr(z[i,genloc[i]] == 1)


# m.addConstr(z['1','12'] ==  1)
# m.addConstr(z['2','32'] == 1)
# m.addConstr(z['3','1'] == 1)
# m.addConstr(z[,78] == 1)
#m.addConstr(z[4,10] == 1)



###########################################################################################
###########################################################################################

print ('test9')


############for microgrid(2)

m.addConstrs(quicksum(v[i,k] for k in K) == 1 for i in BB )
m.addConstrs(v[i,k] >= z[k,i] for i in BB for k in K )

for i,j in list(set(AB) - set(linkoAB)):
	m.addConstrs(c[i,j,k] <= v[i,k] for k in K)
	m.addConstrs(c[i,j,k] <= v[j,k] for k in K)
	m.addConstrs(c[i,j,k] >= v[i,k] + v[j,k] - 1 for k in K)


for i,j in AB:
	m.addConstr(b[i,j] == quicksum(c[i,j,k] for k in  K))


#########################################

###############################For dispatchable 

m.addConstrs(dd[int(i) -1 ]*hp*p[int(i) -1] + (1 - dd[int(i) -1])*p[int(i)-1] <= pd[i] for i in BB )

m.addConstrs( pd[i] <= p[int(i)-1] for i in  BB )

m.addConstrs(ps[i] <= s[i]*p[int(i)-1] for i in BB )

m.addConstrs(ps[i] <= pd[i] for i in BB )

m.addConstrs(ps[i] >= pd[i] - (1- s[i])*p[int(i)-1] for i in BB )

for i in BB:
	if p[int(i)-1] != 0:
		m.addConstr(qs[i] == q[int(i)-1]/p[int(i)-1]*ps[i] )

print ('test10')


########################################################################### disruption
m.addConstrs(s[i] == 0 for i in nodeo)

m.addConstrs(s[i] == 1 for i in nodec)


for item in linkoAB:
	m.addConstr(b[item] == 0)


for item in linkcAB:
	m.addConstr(b[item] == 1)


###########################################################################################
###########################################################################################




print ('test11')

############     OBJECTIVE VALUE    ##############
m.setObjective( quicksum(w[int(i)-1]*ps[i] for i in BB) , GRB.MAXIMIZE)

m.Params.IntFeasTol = 1e-05
m.Params.FeasibilityTol = 1e-06
m.Params.OptimalityTol = 1e-06
m.Params.MIPGap = 1e-02

################################################

m.optimize()
print('\nCost: %g' % m.objVal)


served = 0

for i in m.getVars():
	# print('%s %g' % (i.varName, i.x))
	if i.varName.split('[')[0] == 'location':
		if int(i.x) == 1:
			print('%s %g' % (i.varName, i.x))
	# if i.varName.split('[')[0] == 'activepowergen':
	# 	# if int(i.x) > 0:
	# 	print('%s %g' % (i.varName, i.x))
	if i.varName.split('[')[0] == 'activepowergen':
		if int(i.x) > 0:
			print('%s %g' % (i.varName, i.x))
			# if i.varName.split(',')[1].split(']')[0] == str(j):
			served = i.x + served
			# print('%s %g' % (i.varName, i.x))
	# if i.varName.split('[')[0] == 'reactivepowergen':
	# 	if int(i.x) > 0:
	# 		print('%s %g' % (i.varName, i.x))
	# if i.varName.split('[')[0] == 'activepower':
	# 	if int(i.x) != 0:
	# 		print('%s %g' % (i.varName, i.x))


# all_served = 0
# for i in range(S):
# 	all_served = all_served + served[i]*prob[i]
rae = (served/sum(Pmax))*100

print ('used generation:')
print (served)
print ('RAE(%):')
print (rae)




