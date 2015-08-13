#!/usr/bin/env python2


setpoints = []

values = []

POINTS = False
FIELD = False

with open("postProcessing/sets/modulwinkel/10000/plane_T.vtk", "rb") as f:
    for row in f:
        row.strip()
        if row.startswith("POINTS"):
            print "p"
            POINTS = True
            
        elif row.startswith("POINT_DATA"):
            print "-p"
            POINTS = False
        
        elif row.startswith("T"):
            print "t"
            FIELD = True
        
        elif POINTS:
            setpoints.append([round(float(x), 6) for x in row.strip().split()])
#            setpoints.append(row)
        
        elif FIELD:
            print "yes"
            values = [round(float(x), 6) for x in row.strip().split()]
        
print setpoints
print values

POINTS = False
FIELD = False

surfpoints = []

with open("postProcessing/surfaces/10000/T_plane.vtk", "rb") as f:
    for row in f:
        if row.startswith("POINTS"):
            POINTS = True
            
        elif row.startswith("POLYGONS") or not row.strip():
            POINTS = False
            break
        
        elif POINTS:
            surfpoints.append([round(float(x), 6) for x in row.strip().split()])
#            surfpoints.append(row)

print surfpoints

    
for point in surfpoints:
    if point in setpoints:
        print values[setpoints.index(point)]


