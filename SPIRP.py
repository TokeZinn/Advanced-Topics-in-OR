# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 10:14:08 2019

@author: tokez
"""

#Import Required Packages
import math
import numpy
import matplotlib.pyplot as plt 
from gurobipy import *
from scipy.stats import binom
import time
from matplotlib.pyplot import cm
from tqdm import tqdm
from copy import copy

#Full runtime calculation
LARGE_T = time.time()

#Initialize data
instance = numpy.loadtxt(fname = "data_SPIRP.txt")


## ---------- GENERATION VARIABLES ---------- ## 

# -- Simulation Setup -- #
MaxIter = 1 # Number of Simulations
Periods = 1  # T = 30
nStores = 40 # S = 40

# -- Model Parameters -- #
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

# -- Print Settings -- #
Draw = True # See graph/plot of solution for each time
PrintSummary = False #If a summary should be printed

# -- Our Additions -- # 
AllowPartial = True #Allow a node to be visited on multiple routes. 
MyopicPolishing = True #Apply the myopic polishing procedure - Algorithm 4



## ---------------- FUNCTIONS --------------- ##

#Distance between nodes.
def Distance(i,j):
    return(math.sqrt((instance[i][1]-instance[j][1])**2 + 
                     (instance[i][2]-instance[j][2])**2))
    
# Obtain delivery quantity based on inventory
def Delivery(Inventory):
    if Inventory >= EV:
        return(0)
    else:
        return(min(C-Inventory, numpy.floor((L)*EV) - Inventory))

#The final point of a route/path            
def End(path):
    return(path[len(path)-1])

#Route Extension Function
def ExtendPath(path,maxpath = 10000):
    local = ()
    #For each store that needs delivery
    for i in active:
        #if the store could not be added to subpath, is in path or is appended continue
        if i in path["n"] or i in path["p"] or i in path["a"]: 
            continue
        else: 
            D = path["d"] + Dist[End(path["p"]),i] + Dist[i,0]
            #if the distance of the route remains feasible and we have sufficient capcity
            if D <= maxTrip and path["Qty"] + quantities[i] <= Q + AllowPartial*(-1+min(x for x in quantities.values() if x>0)):
                #add the store to be added
                local = local + (i,)
            else:
                #Else save the store so that paths generated from this does not need to check 
                # if the node is compatible
                path["n"] = path["n"] + (i,)
    for i in local:
        #Add the found notes if they have not already been added
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

#Update the Inventory
def UpdateInventory(inventory, demand, quantities):
    #Coerce to list
    delivery = []
    for i in quantities:
        delivery += [quantities[i],]
    
    #Add the delivery to the "freshest" row    
    inventory[L-1,:] += delivery
    
    #Subtract the demand
    for s in N:
        #first from the "oldest" inventory and then cascade to new inventory
        for l in range(L):    
            if inventory[l,s] >= demand[s]:
                inventory[l,s] = inventory[l,s] - demand[s]
                demand[s] = 0
            else:
                demand[s] = demand[s] - inventory[l,s]
                inventory[l,s] = 0
    #Make a copy of the current row as it is waste. 
    #NOTE: copy() is important here as a pointer is otherwise stored and not 
    #called til the return command is executed. Hence resulting in the wrong
    #amount of waste!         
    waste = copy(inventory[0,:])
    inventory[0,:] = numpy.zeros((1,len(inventory[0,:])))
    inventory[tuple(range(L)),:] = inventory[tuple(range(1,L)) + (0,),:]
    return(waste)
    
## ------------ MYOPIC PROCEDURES ----------- ##

# Immediate Profit Calculation
def ImmediateProfit(i,route):
    c_r = sum(Dist[route[j],route[j+1]] for j in range(len(route)-1))
    g = tuple(x for x in route if x != i)
    c_g = sum(Dist[g[j],g[j+1]] for j in range(len(g)-1))
    delta = c_r - c_g 
    
    ImProf = (p-a)*((EV) - initial[i]) - delta
    return(ImProf)
    



## ------------ START GENERATION ------------ ##
    
# Initialize graph
Dist = {(i,j): Distance(i,j) for i in N for j in N}
Results = {}
Waste = {}


# Capacity of Store varies with Shelf-Life
if L == 2:
    initialrange = 30+1
    C = 40
if L == 3:
    initialrange = 50+1
    C = 60
if L == 4:
    initialrange = 70+1
    C = 80


#Set seed for reproducability
numpy.random.seed(101)

#Used for storing time used in each run
SIM_T = 0

for iteration in tqdm(range(MaxIter)):
    if not AllowPartial and MyopicPolishing:
        print("Runnng Myopic Polishing without Partial Deliveries is NOT recommended.")
    #iteration timer
    tbar = time.time()
    
    #Initilize Simulation Dictionaries for storing results. 
    VRP = {}
    R = {}
    Z = {}
    Y = {}
    D = {}
    Demand = {}
    Inventory = {}
    Solution = {}
    profit = {}
    waste = {}
    aquired = {}
    
    #Initialize inventory (note that the depot is also initialized but set to 0)
    initial = [numpy.random.randint(0,initialrange) for i in N]
    initial[0] = 0
    first = sum(j for j in initial)
    inventory = numpy.zeros((L,nStores+1))
    inventory[L-2,:] = initial
    
    for tau in range(1,Periods+1):
        if Draw:
            print('--------------------\\\\\--------------------')
            print('                  Period', tau)
            print('--------------------\\\\\--------------------')
            print()
        
        #Initialize trucks, demand etc. Note demand[0] is set to 0, as 0 is the depot
        Trucks = {}
        demand = numpy.random.binomial(200, 0.1, size=nStores+1)
        demand[0] = 0
        Demand[tau] = demand
        Inventory[tau] = inventory
        
        
        #First we decide whether or not we want to deliver 
        active = []
        for i in range(1,nStores): 
            if initial[i] < EV:
                active += [i,]
        
        quantities = {(i): Delivery(initial[i]) for i in N}
        quantities[0] = 0
        
        #Initialize trucks
        Trucks = {i: {"p": (0,active[i]), "d": Dist[0,active[i]], "Qty": quantities[active[i]], "n": (),"a": ()} 
            for i in range(len(active)) if Dist[0,active[i]] <= maxTrip}
        
        #Generate set of trucks for t: R_t
        t = time.time()
        start = 0
        while True:
            Len = len(Trucks)
            for i in range(start,len(Trucks)):
                ExtendPath(Trucks[i])
            if Len == len(Trucks):
                break
            else:
                start = Len 
        
        #If we are drawing 
        if Draw:
            print('Generating trucks took: %s' % (time.time() - t))
            print()
        
        #local timer
        t = time.time()
        
        #Initilize variables
        VRP[tau] = Model('VRP')
        VRP[tau].setParam('OutputFlag',0)
        R[tau] = {r: Trucks[r]["p"] + (0,) for r in range(len(Trucks))}
        Z[tau] = {r: VRP[tau].addVar(vtype = GRB.BINARY) for r in R[tau]}
        Y[tau] = {(r,i): VRP[tau].addVar(vtype= GRB.INTEGER) for r in R[tau] for i in active if i in R[tau][r]}
        D[tau] = {r: Trucks[r]["d"] + Dist[End(Trucks[r]["p"]),0] for r in R[tau]}
        
        #Setup problem - see project for description
        VRP[tau].setObjective(quicksum(Z[tau][r]*D[tau][r] for r in R[tau]), GRB.MINIMIZE)
        DeliveryIfRoute = {r: VRP[tau].addConstr(quicksum(Y[tau][r,i] for i in R[tau][r] if i > 0) <= Q*Z[tau][r]) for r in R[tau]}
        DeliverQuantity = {i: VRP[tau].addConstr(quicksum(Y[tau][r,i] for r in R[tau] if i in R[tau][r]) == quantities[i]) for i in active}
        
        #Set branch priorit of the routes. 
        for r in R[tau]:
            Z[tau][r].BranchPriority = 1
        
        #Optimize the problem and collect the routes that are used
        VRP[tau].optimize()
        R_star = {r: R[tau][r] for r in Z[tau] if Z[tau][r].x > 0.5}
        
        
        #Draw Plots
        if Draw and not MyopicPolishing:
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
                   plt.plot(x,y, marker = "o",linestyle = ":" if MyopicPolishing else "-",c = color[j-1,:] )
            
            plt.plot(instance[active,1],instance[active,2],"ro", color = "green")
            
            plt.plot(instance[0,1],instance[0,2], "ro", color = "red")
            if not MyopicPolishing:
                plt.show()
        
        # Apply the Myopic Polishing Procedure
        if MyopicPolishing:
            #For each r and node in r, test if it can be and if it reduces
            #the cost more than the current node. Break if no node is found.
            #Then overwrite the current solution. 
            for r in R_star:
                while True:
                    i_min = None
                    minimum = 0
                    for i in R_star[r]:
                        if i == 0:
                            continue
                        else:
                            if ImmediateProfit(i,R_star[r]) < minimum:
                                i_min = i
                                minimum = ImmediateProfit(i,R_star[r])
                    if (i_min == None):
                        break
                    else:
                        R_star[r] = tuple(x for x in R_star[r] if x != i_min)
                        quantities[i_min] = quantities[i_min] - Y[tau][r,i_min].x
            R_star = {r: R_star[r] for r in R_star if R_star[r] != (0,0)}
        
        
        # Graph plot of Myopic Solution                
        if Draw & MyopicPolishing:
            color=cm.rainbow(numpy.linspace(0,1,len(R_star)))
            plt.plot(instance[range(nStores),1],instance[range(nStores),2],"ro", color = "grey")
            plt.plot(instance[active,1],instance[active,2],"ro", color = "green")

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
            
                            
        
        
        #Calculate Sold Units
        sold_units = sum(min(demand[i],initial[i] + quantities[i]) for i in range(1,nStores+1))
        
        # Calculate Acquired Units
        aquired_units = sum(quantities[i] for i in range(1,nStores+1))
        aquired[tau] = aquired_units
        
        #Calculate Waste
        temp = UpdateInventory(inventory,demand,quantities)
        waste[tau] = sum(i for i in temp)
        
        #Calculate Aggregated Inventory
        initial = numpy.sum(inventory,0)
        
        #Calculate Distance Cost of Solution
        distance_cost = sum(sum(Dist[R_star[r][i],R_star[r][i+1]] for i in range(len(R_star[r])-1)) for r in R_star)
        
        #Calculate Profit
        profit[tau] = p*sold_units - a*aquired_units - distance_cost
        
    #calculate time cost of iteration
    SIM_T = SIM_T + time.time()-tbar
    
    if PrintSummary:
        print('--------------------\\\\\--------------------')
        print('           Summary for Iteration %d' %(iteration))
        print('--------------------\\\\\--------------------')
        print()
        
        
        
        print("Total time elapsed: %s seconds" %(time.time()-LARGE_T))
        print()
        print("Total profit: %s" %(sum(profit[tau] for tau in range(1,Periods+1))))
        print()
    
    #Save iteration results
    Results[iteration] = sum(profit[tau] for tau in range(1,Periods+1))
    total = sum(aquired[tau] for tau in range(1,Periods+1)) + first
    Waste[iteration] = sum(waste[tau] for tau in range(1,Periods+1))/total
    
#Print result of simulations.
print()
print('--------------------\\\\\--------------------')
print('        Simulation Results for L = %s' %(L))
print('--------------------\\\\\--------------------')
print()    
if MyopicPolishing:
    print("Solution was polished to be more myopic")
if AllowPartial:
    print("Multiple visits to a store in different routes allowed")
print("Average Profit over %s runs: %d" %(MaxIter,numpy.mean(list(Results.values()))))
print("Average Waste over %s runs: %.5f" %(MaxIter,numpy.mean(list(Waste.values()))))
print("Average Runtime: %.6f" %(SIM_T/MaxIter))

