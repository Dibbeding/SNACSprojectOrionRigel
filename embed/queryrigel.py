#!/usr/bin/env python

import math
import sys

def query(coordinates, source, destination, curvature):
    def arccosh(x):
        t = x + math.sqrt(x * x - 1)
        return math.log(t)

    sourceCoords = coordinates[source]
    destinationCoords = coordinates[destination]

    i = 0
    ts = 1.0
    td = 1.0
    tt = 1.0

    for i in range(len(sourceCoords)):
        ts += math.pow(sourceCoords[i], 2)
        td += math.pow(destinationCoords[i], 2)
        tt += (sourceCoords[i] * destinationCoords[i])

    t = math.sqrt(ts * td) - tt
    return arccosh(t) * math.fabs(curvature)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit('Usage: python %s [coordinate file prefix] [# of partitions] [curvature]' % sys.argv[0])

    landFile = sys.argv[1] + '.land'
    partitions = int(sys.argv[2])
    coordFiles = [sys.argv[1] + str(i) + '.coord' for i in range(partitions)]
    curvature = int(sys.argv[3])

    coordinates = dict()
    with open(landFile) as infile:
        for line in infile:
            linesplit = line.split()
            id = int(linesplit[0])
            coords = [float(c) for c in linesplit[1:]]
            coordinates[id] = coords

    for coordFile in coordFiles:
        with open(coordFile) as infile:
            for line in infile:
                linesplit = line.split()
                id = int(linesplit[0])
                coords = [float(c) for c in linesplit[1:]]
                coordinates[id] = coords

    while True:
        query_input = raw_input("Enter ID of 2 nodes: ")
        
        if query_input == 'exit' || query_input == 'q' || query_input == 'quit':
            break

        querysplit = query_input.split()
        source = int(querysplit[0])
        destination = int(querysplit[1])

        estimate = query(coordinates, source, destination, curvature)
        print 'Rigel estimates the distance between %d and %d to be %f.\n' % (source, destination, estimate)
