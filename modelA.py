
import sys
import string

from ortools.constraint_solver import pywrapcp

# This function create the constraints that links the position of the lower left vertex to the one of the upper right vertex
# Parameters:
#
# solver:         the constraint solver
# xOrigin:        array of IntVar
# xDestination:   array of IntVar
# yOrigin:        array of IntVar
# yDestination:   array of IntVar
# zOrigin:        array of IntVar
# zDestination:   array of IntVar
# length:         array of int
# width:          array of int
# height:         array of int
def correspondingOriginDestination(solver,xOrigin,xDestination,yOrigin,yDestination,zOrigin,zDestination,length,width,height):
  for i in range(len(length)):
    if length[i] == width[i]:
      solver.Add(xOrigin[i] + length[i] == xDestination[i])
      solver.Add(yOrigin[i] + length[i] == yDestination[i])
    else:
      c1 = xOrigin[i] + length[i] == xDestination[i]
      c2 = xOrigin[i] + width[i] == xDestination[i]
      solver.Add(c1+c2 == 1)
      c3 = yOrigin[i] + length[i] == yDestination[i]
      c4 = yOrigin[i] + width[i] == yDestination[i]
      solver.Add(c3+c4 == 1)
      solver.Add(yDestination[i]-yOrigin[i]+xDestination[i]-xOrigin[i] == length[i]+width[i])

    solver.Add(zOrigin[i] + height[i] == zDestination[i])

# This function create the constraints that forbid the box from intersecting.
# Parameters:
#
# solver:         the constraint solver
# xOrigin:        array of IntVar
# xDestination:   array of IntVar
# yOrigin:        array of IntVar
# yDestination:   array of IntVar
# zOrigin:        array of IntVar
# zDestination:   array of IntVar
def noIntersection(solver,xOrigin,xDestination,yOrigin,yDestination,zOrigin,zDestination):
  for i in range(len(xOrigin)):
    for j in range(i+1,len(xOrigin)):
      c1 = xOrigin[i] >= xDestination[j]
      c2 = xOrigin[j] >= xDestination[i]
      c3 = yOrigin[i] >= yDestination[j]
      c4 = yOrigin[j] >= yDestination[i]
      c5 = zOrigin[i] >= zDestination[j]
      c6 = zOrigin[j] >= zDestination[i]
      solver.Add(c1+c2+c3+c4+c5+c6 >= 1)

# This function create the constraints that force the box to be either on the ground or fully on a single other box.
# Parameters:
#
# solver:         the constraint solver
# xOrigin:        array of IntVar
# xDestination:   array of IntVar
# yOrigin:        array of IntVar
# yDestination:   array of IntVar
# zOrigin:        array of IntVar
# zDestination:   array of IntVar
def onSurface(solver, xOrigin, xDestination, yOrigin, yDestination, zOrigin, zDestination):
  for i in range(len(xOrigin)):
    bb = []
    c1 = zOrigin[i] == 0
    bb.append(c1)
    for j in range(len(xOrigin)):
      if not i == j:
        c2 = xOrigin[i] >= xOrigin[j]
        c3 = xDestination[i] <= xDestination[j]
        c4 = yOrigin[i] >= yOrigin[j]
        c5 = yDestination[i] <= yDestination[j]
        c6 = zOrigin[i] == zDestination[j]
        bb.append(c2*c3*c4*c5*c6)
    solver.Add(solver.Sum(bb) == 1)


# This function create the constraints that indicate in which stop the box are unloaded, when a box is unloaded it has to be accessible.
# Parameters:
#
# solver:         the constraint solver
# xOrigin:        array of IntVar
# xDestination:   array of IntVar
# yOrigin:        array of IntVar
# yDestination:   array of IntVar
# zOrigin:        array of IntVar
# zDestination:   array of IntVar
# groups:         array of integer indicating in which stop the box has to be unloaded
def forceOrder(solver, xOrigin, xDestination, yOrigin, yDestination, zOrigin, zDestination, groups):
  for i in range(len(xOrigin)):
    for j in range(len(xOrigin)):
      if groups[j] < groups[i]:
        ci1 = xOrigin[j] >= xDestination[i] # Sufficient condition
        ci2 = yOrigin[j] >= yDestination[i] # Sufficient condition
        ci3 = yOrigin[i] >= yDestination[j] # Sufficient condition
        ci4 = xOrigin[j] >= xOrigin[i]
        ci5 = zOrigin[j] >= zDestination[i] # ci4 && ci5 Sufficient condition
        solver.Add(ci1 + ci2 + ci3 + (ci4*ci5) >= 1)

def breakSymmetries(solver, xOrigin, xDestination, yOrigin, yDestination, zOrigin, zDestination, groups, length, width, height):
  for i in range(len(length)):
    for j in range(i+1,len(length)):
      if groups[i] == groups[j] and length[i] == length[j] and width[i] == width[j] and height[i] == height[j]:
        ci1 = xOrigin[i] < xOrigin[j]
        ci2 = xOrigin[i] == xOrigin[j]
        ci3 = yOrigin[i] < yOrigin[j]
        ci4 = yOrigin[i] == yOrigin[j]
        ci5 = zOrigin[i] < zOrigin[j]
        solver.Add(ci1 + (ci2*ci3) + (ci2*ci4*ci5) >= 1) 

# Parse the instance
#
# Parameters:
# path: the path to the instance file.
#
# Return:
# the data of the instance (container dimensions, boxes)
def parseInstance(path):
  f = open(path,'r')
  lines = f.readlines()
  container = [int(lines[0].split(' ')[0]),int(lines[0].split(' ')[1]),int(lines[0].split(' ')[2])]
  boxes = []
  for line in lines[1:]:
    boxes.append([int(line.split(' ')[0]),int(line.split(' ')[1]),int(line.split(' ')[2]),int(line.split(' ')[3])])
  return (container, boxes)

# Build and solve the model
#
# Parameters:
#
# instance the data of the problem (container dimensions, boxes)
def buildAndSolveModel(instance):
  solver = pywrapcp.Solver('Container Loader')
  container = instance[0]
  boxes = instance[1]
  boxesInst = [[box[i] for box in boxes] for i in range(4)]
  length = boxesInst[0]
  width = boxesInst[1]
  height = boxesInst[2]
  order = boxesInst[3]
  xOrigin = [solver.IntVar(0,container[0],'xOrigin[%i]'%i) for i in range(len(length))]
  xDestination = [solver.IntVar(0,container[0],'xDestination[%i]'%i) for i in range(len(length))]
  yOrigin = [solver.IntVar(0,container[1],'yOrigin[%i]'%i) for i in range(len(length))]
  yDestination = [solver.IntVar(0,container[1],'yDestination[%i]'%i) for i in range(len(length))]
  zOrigin = [solver.IntVar(0,container[2],'zOrigin[%i]'%i) for i in range(len(length))]
  zDestination = [solver.IntVar(0,container[2],'zDestination[%i]'%i) for i in range(len(length))]
  correspondingOriginDestination(solver,xOrigin,xDestination,yOrigin,yDestination,zOrigin,zDestination,length,width,height)
  noIntersection(solver,xOrigin,xDestination,yOrigin,yDestination,zOrigin,zDestination)  
  onSurface(solver,xOrigin,xDestination,yOrigin,yDestination,zOrigin,zDestination)
  forceOrder(solver, xOrigin, xDestination, yOrigin, yDestination, zOrigin, zDestination, order)
  breakSymmetries(solver, xOrigin, xDestination, yOrigin, yDestination, zOrigin, zDestination, order, length, width, height)


  solution = solver.Assignment()
  solution.Add(xOrigin)
  solution.Add(xDestination)
  solution.Add(yOrigin)
  solution.Add(yDestination)
  solution.Add(zOrigin)
  solution.Add(zDestination)

  var = []
  orderedGroups = [(order[i],length[i]*width[i],i) for i in range(len(length))]
  orderedGroups.sort(reverse=True)
  for (gp,dummy,i) in orderedGroups:
    var.append(xOrigin[i])
    var.append(xDestination[i])
    var.append(yOrigin[i])
    var.append(yDestination[i])
    var.append(zOrigin[i])
    var.append(zDestination[i])
  objVar =  solver.Max(xDestination)
  objective = solver.Minimize(objVar,1)

  db = solver.Phase(var, solver.CHOOSE_FIRST_UNBOUND, solver.ASSIGN_MIN_VALUE)
  solver.NewSearch(db, [objective])

  print "Start Search"
  while solver.NextSolution():
    print "xOrigin: ", [xOrigin[i].Value() for i in range(len(xOrigin))]
    print "yOrigin: ", [yOrigin[i].Value() for i in range(len(xOrigin))]
    print "zOrigin: ", [zOrigin[i].Value() for i in range(len(xOrigin))]
    print "xDestination: ", [xDestination[i].Value() for i in range(len(xOrigin))]
    print "yDestination: ", [yDestination[i].Value() for i in range(len(xOrigin))]
    print "zDestination: ", [zDestination[i].Value() for i in range(len(xOrigin))]
    print "----------"
  solver.EndSearch()
  print "End Search"
  print "Failures:", solver.Failures()
  print "Branches:", solver.Branches()
  print "WallTime:", solver.WallTime()

# Main function
# Expects arguments:
# 1/ Path to the instance file
def main():
  args = sys.argv[1:]
  instance = parseInstance(args[0])
  
  buildAndSolveModel(instance)

if __name__ == '__main__':
  main()
