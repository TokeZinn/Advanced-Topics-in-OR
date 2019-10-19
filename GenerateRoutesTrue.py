
import math
import numpy
import matplotlib.pyplot as plt 
from gurobipy import *
from scipy.stats import binom
import time
from matplotlib.pyplot import cm
from tqdm import tqdm




LARGE_T = time.time()

instance = numpy.loadtxt(fname = "data_SPIRP.text")

def FullProbSize(n):
    return(sum(math.factorial(i) for i in range(2,n+1)))

def Distance(i,j):
    return(math.sqrt((instance[i][1]-instance[j][1])**2 + 
                     (instance[i][2]-instance[j][2])**2))
    
def Delivery(Inventory):
    if Inventory >= EV:
        return(0)
    else:
        return(min(C-Inventory, numpy.floor((L)*EV) - Inventory))
            
def End(path):
    return(path[len(path)-1])
        
def ExtendPath(path, maxpath = 10000, heuristics = False):
    local = ()
    for i in active:
        if i in path["n"] or i in path["p"] or i in path["a"]: 
            continue
        else: 
            D = path["d"] + Dist[End(path["p"]),i] + Dist[i,0]
            if heuristics: 
                if Dist[End(path["p"]),i] >= MaxSingleDist:
                    continue
            if D <= maxTrip and path["Qty"] + quantities[i] <= Q + AllowPartial*(-1+min(x for x in quantities.values() if x>0)):
                local = local + (i,)
            else:
                path["n"] = path["n"] + (i,)
    for i in local:
        if(i in path["a"]):
            continue
        else:
            if len(path["p"] + (i,)) <= maxpath:
                Trucks[len(Trucks)] = {"p": path["p"] + (i,),
                       "d": path["d"] + Dist[End(path["p"]),i],
                       "Qty": path["Qty"] + quantities[i],
                       "n": path["n"], "a": ()}
                path["a"] += (i,)
            else:
                continue

def UpdateInventory(inventory, demand, quantities):
    delivery = []
    for i in quantities:
        delivery += [quantities[i],]
        
    inventory[L-1,:] += delivery
    
    for l in range(L):
        for s in N:
            inventory[l,s] = max(0, inventory[l,s] - demand[s])
            demand[s] = max(0, demand[s] - inventory[l,s])
    
    inventory[0,:] = numpy.zeros((1,len(inventory[0,:])))
    inventory[tuple(range(L)),:] = inventory[tuple(range(1,L)) + (0,),:]
    


## ------------ MYOPIC PROCEDURES ----------- ##
    
def ImmediateProfit(i,route,shift = 0):
    c_r = sum(Dist[route[j],route[j+1]] for j in range(len(route)-1))
    g = tuple(x for x in route if x != i)
    c_g = sum(Dist[g[j],g[j+1]] for j in range(len(g)-1))
    delta = c_r - c_g 
    
    ImProf = (p-a)*((EV+shift) - initial[i]) - delta
    return(ImProf)
    
    
def IntersectingRoutes(R):
    H = {}
    for r in R:
        for r_dash in R:
            if r == r_dash:
                continue
            else:
                for i in active:
                    if i in R[r] and i in R[r_dash]:
                        H[i] = (r,r_dash)
    return(H)


## ---------- LONG-TERM PROCEDURES ---------- ##
    
def DelayDelivery(i,Q_b):
    return(((initial[i] + Q_b) // EV) - (initial[i] // EV))
    
def FindLeastForDelay(i,Q_b):
    for qty in range(1,Q_b+1):
        if DelayDelivery(i,qty) == DelayDelivery(i,Q_b):
            return(qty)
    
def LongTermProfit(i,route,qty, shift = 0):
    c_r = sum(Dist[route[j],route[j+1]] for j in range(len(route)-1))
    c_rp = math.inf
    r_best = None
    for j in range(1,len(route)-1):
        r = route
        r = r[:j] + (i,) + r[j:]
        c_test = sum(Dist[r[k],r[k+1]] for k in range(len(r)-1))
        if c_test < c_rp:
            c_rp = c_test
            r_best = r
    delta = c_rp - c_r
    return({i: {"r": r_best, "profit": (p-a)*qty - delta }})
    
        


## ---------- GENERATION VARIABLES ---------- ##

MaxPath = 1000
LimitPath = True
MaxSingleDist = 400
Heuristics = False
Draw = False
AllowPartial = False
PrintSummary = False
MaxIter = 30
Periods = 30  # T = 30
nStores = 40 # S = 40
T = range(Periods+1) #Periods
N = range(nStores+1) #Number of Vertices
S = range(1,nStores+1) #Number of store-only vertices
L = 4 # deterministic shelf life in {2,3,4}
a = 6 # acquisition cost
p = 10 # selling price
EV = 20 # expected value of demand
Q = 120 # vehicle capacity
maxTrip = 230 # maximum trip length
TSL = 0.9 # target service level

LongTermPolishing = False

MyopicPolishing = False
level = 0.5 #... 
P_shift = binom.ppf(level, 200,0.1) - EV


Dist = {(i,j): Distance(i,j) for i in N for j in N}
Arcs = {(i,j) for i in N for j in N if j != i}

## ------------ START GENERATION ------------ ##


Results = {}

if L == 2:
    initialrange = 30+1
    C = 40
if L == 3:
    initialrange = 50+1
    C = 60
if L == 4:
    initialrange = 70+1
    C = 80

numpy.random.seed(100)

SIM_T = 0

for iteration in tqdm(range(MaxIter)):
    VRP = {}
    R = {}
    Z = {}
    Y = {}
    P = {}
    D = {}
    Demand = {}
    Inventory = {}
    Solution = {}
    profit = {}
    
    initial = [numpy.random.randint(0,initialrange) for i in N]
    initial[0] = 0
    #demand = [numpy.random.binomial(200, 0.1, size=nStores+1) for t in T]
    inventory = numpy.zeros((L,nStores+1))
    inventory[L-2,:] = initial
    
    for tau in tqdm(range(1,Periods+1)):
        tbar = time.time()
        if Draw:
            print('--------------------\\\\\--------------------')
            print('                  Period', tau)
            print('--------------------\\\\\--------------------')
            print()
        Trucks = {}
        demand = numpy.random.binomial(200, 0.1, size=nStores+1)
        Demand[tau] = demand
        Inventory[tau] = inventory
        
        
        #First we decide whether or not we want to deliver 
        active = []
        for i in range(1,nStores): 
            if initial[i] < EV:
                active += [i,]
        
        quantities = {(i): Delivery(initial[i]) for i in N}
        quantities[0] = 0
        
    
        Trucks = {i: {"p": (0,active[i]), "d": Dist[0,active[i]], "Qty": quantities[active[i]], "n": (),"a": ()} 
            for i in range(len(active)) if Dist[0,active[i]] <= maxTrip}
        
        t = time.time()
        start = 0
        while True:
            Len = len(Trucks)
            for i in range(start,len(Trucks)):
                ExtendPath(Trucks[i],MaxPath*LimitPath,Heuristics)
            if Len == len(Trucks):
                break
            else:
                start = Len 
        
        if Draw:
            print('Generating trucks took: %s' % (time.time() - t))
            print()
        
        t = time.time()
        VRP[tau] = Model('VRP')
        VRP[tau].setParam('OutputFlag',0)
        R[tau] = {r: Trucks[r]["p"] + (0,) for r in range(len(Trucks))}
        Z[tau] = {r: VRP[tau].addVar(vtype = GRB.BINARY) for r in R[tau]}
        Y[tau] = {(r,i): VRP[tau].addVar(vtype= GRB.INTEGER) for r in R[tau] for i in active if i in R[tau][r]}
        P[tau] = {r: Trucks[r]["Qty"]*(p-a) for r in R[tau]}
        D[tau] = {r: Trucks[r]["d"] + Dist[End(Trucks[r]["p"]),0] for r in R[tau]}
        
        VRP[tau].setObjective(quicksum(Z[tau][r]*D[tau][r] for r in R[tau]), GRB.MINIMIZE)
        DeliveryIfRoute = {r: VRP[tau].addConstr(quicksum(Y[tau][r,i] for i in R[tau][r] if i > 0) <= Q*Z[tau][r]) for r in R[tau]}
        DeliverQuantity = {i: VRP[tau].addConstr(quicksum(Y[tau][r,i] for r in R[tau] if i in R[tau][r]) == quantities[i]) for i in active}
        
        for r in R[tau]:
            Z[tau][r].BranchPriority = 1
        
        VRP[tau].optimize()
        R_star = {r: R[tau][r] for r in Z[tau] if Z[tau][r].x > 0.5}
        
        if Draw:
            print('Solving VRP took: %s' % (time.time() - t))
            print()
            color=cm.rainbow(numpy.linspace(0,1,len(R_star)))
            plt.plot(instance[range(nStores),1],instance[range(nStores),2],"ro", color = "grey")
            plt.plot(instance[active,1],instance[active,2],"ro", color = "green")
            plt.plot(instance[0,1],instance[0,2], "ro", color = "red")
            
            j = 0
            for r in R_star:
                #Add each route in the solution
                j+= 1
                q = R_star[r]
                for i in range(len(q)-1):
                    x = [instance[q[i],1], instance[q[i+1],1]]
                    y = [instance[q[i],2], instance[q[i+1],2]]
                    plt.plot(x,y, marker = "o",linestyle = ":",c = color[j-1,:] )
            
            plt.plot(instance[active,1],instance[active,2],"ro", color = "green")
            
            plt.plot(instance[0,1],instance[0,2], "ro", color = "red")
            if not MyopicPolishing:
                plt.show()
        
        
        if MyopicPolishing:
            print("Entered")
            for r in R_star:
                while True:
                    i_min = None
                    minimum = 0
                    for i in R_star[r]:
                        if i == 0:
                            continue
                        else:
                            if ImmediateProfit(i,R_star[r],P_shift) < minimum:
                                i_min = i
                                minimum = ImmediateProfit(i,R_star[r],P_shift)
                    if (i_min == None):
                        break
                    else:
                        R_star[r] = tuple(x for x in R_star[r] if x != i_min)
                        quantities[i_min] = quantities[i_min] - Y[tau][r,i_min].x
            R_star = {r: R_star[r] for r in R_star if R_star[r] != (0,0)}
        

        if LongTermPolishing:
            Q_bar = {r: Q - sum(Y[tau][r,i].x for i in R_star[r] if i in active) for r in R_star}
            Qty = {}
            for r in Q_bar:
                if Q_bar[r] < 0.1:
                    continue
                else:
                    for i in range(1,nStores+1):
                        if (i in R_star[r]) or (i in active) or initial[i] == C:
                            continue
                        else:
                            if (DelayDelivery(i,Q_bar[r]) > 0):
                                q_ri = FindLeastForDelay(i,int(Q_bar[r]))
                                if initial[i] + q_ri <= C:
                                    Qty[r,i] = q_ri                    
            for r in Q_bar:
                Q_bar[r] = int(Q_bar[r])
            while True:
                p_best = 0
                r_best = None
                i_best = None 
                route_best = None
                for (r,i) in Qty:
                    temporary = LongTermProfit(i,R_star[r],Qty[r,i])
                    if (temporary[i]["profit"] > p_best) and (Q_bar[r] - Qty[r,i] >= 0):
                        p_best = temporary[i]["profit"]
                        r_best = r
                        i_best = i
                        route_best = temporary[i]["r"]
                if i_best == None:
                    break
                
                Q_bar[r_best] += -Qty[r_best,i_best]
                R_star[r_best] = route_best
                initial[i_best] += Qty[r_best,i_best]
                for k in list(Qty.keys()):
                    if k[1] == i_best:
                        Qty.pop(k)
            
                        
        if Draw & MyopicPolishing:
            color=cm.rainbow(numpy.linspace(0,1,len(R_star)))
            plt.plot(instance[range(nStores),1],instance[range(nStores),2],"ro", color = "grey")
            plt.plot(instance[active,1],instance[active,2],"ro", color = "green")
            for i in active:
                plt.annotate(i,(instance[i,1],instance[i,2]))
                plt.plot(instance[0,1],instance[0,2], "ro", color = "red")

            j = 0
            for r in R_star:
                    #Add each route in the solution
                    j+= 1
                    q = R_star[r]
                    for i in range(len(q)-1):
                        x = [instance[q[i],1], instance[q[i+1],1]]
                        y = [instance[q[i],2], instance[q[i+1],2]]
                        plt.plot(x,y, marker = "o",linestyle = "-",c = color[j-1,:] )

            plt.plot(instance[active,1],instance[active,2],"ro", color = "green")

            plt.plot(instance[0,1],instance[0,2], "ro", color = "red")
            plt.show()
            
                            
        
        
        
        sold_units = sum(min(demand[i],initial[i] + quantities[i]) for i in N)
        
        aquired_units = sum(quantities[i] for i in N)
        
        UpdateInventory(inventory,demand,quantities)
        initial = numpy.sum(inventory,0)
        
        #Routes that are used
        distance_cost = sum(sum(Dist[R_star[r][i],R_star[r][i+1]] for i in range(len(R_star[r])-1)) for r in R_star)
        
        profit[tau] = p*sold_units - a*aquired_units - distance_cost
        
        
        #Initialze plot with Data
        
    SIM_T = SIM_T + time.time()-tbar
    potential = a*sum(initial)
    if PrintSummary:
        print('--------------------\\\\\--------------------')
        print('                   Summary')
        print('--------------------\\\\\--------------------')
        print()
        
        
        
        print("Total time elapsed: %s seconds" %(time.time()-LARGE_T))
        print()
        print("Total profit: %s" %(sum(profit[tau] for tau in range(1,Periods+1))))
        print()
        print("Total potential: %s" %(potential + sum(profit[tau] for tau in range(1,Periods+1))))
    Results[iteration] = potential + sum(profit[tau] for tau in range(1,Periods+1))

print()
print('--------------------\\\\\--------------------')
print('        Simulation Results for L = %s' %(L))
print('--------------------\\\\\--------------------')
print()    
if MyopicPolishing:
    print("Solution was polished to be more myopic")
print("Average Profit over %s runs: %d" %(MaxIter,numpy.mean(list(Results.values()))))
print("Average Runtime: %.6f" %(SIM_T/MaxIter))

