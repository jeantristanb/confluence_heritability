#!/usr/bin/python3
import sys, os, gzip,argparse

def parseArguments():
    parser = argparse.ArgumentParser(description='format file for gcta, append N and frequencie if not present using bed file')
    parser.add_argument('--bim',type=str,required=True, help="association files")
    parser.add_argument('--map',type=str,required=True, help="association files")
    parser.add_argument('--chro',type=str,required=True, help="association files")
    parser.add_argument('--out', type=str,help="out of tex file",default="map.map")
    args = parser.parse_args()
    return args

args = parseArguments()


infile = open(args.bim) #.bed
mapfile = gzip.open(args.map) #input map file, either the HapMap map or the 1000 Genomes OMNI map
outfile = gzip.open(args.out, "w") #output style: [rs] [pos] [genetic pos]
chro=args.chro

posin = list()
rsin = list()
mappos = list()
mapgpos = list()

for line in infile :
  line = line.strip().split()
  if line[0] != chro :
     continue
  pos = int(line[3])
  rs = line[1]
  posin.append(pos)
  rsin.append(rs)

infile.close()

line = mapfile.readline()
for line in mapfile:
  line = line.decode('utf-8').strip().split()
  pos = int(line[0])
  gpos = float(line[2])
  mappos.append(pos)
  mapgpos.append(gpos)

mapfile.close()

if len(posin)==0 :
   print('error posin empty')
   sys.exit(2)

def writegz(vec, out):
  vecm=" ".join([str(x) for x in vec])+'\n'
  out.write(vecm.encode())

writeout=outfile
index1 = 0
index2 = 0
while index1 < len(posin):
  pos = posin[index1]
  rs = rsin[index1]
  if pos == mappos[index2]:
    writegz([rs, pos, mapgpos[index2]], writeout)
    index1 = index1+1
  elif pos < mappos[index2]:
    if index2 ==0:
      #before the first site in the map (genetic position = 0)
      
      writegz([rs, pos, mapgpos[index2]], writeout)
      index1 = index1+1
    else:
      #interpolate
      prevg = mapgpos[index2-1]
      prevpos = mappos[index2]
      frac = (float(pos)-float(mappos[index2-1]))/ (float(mappos[index2]) - float(mappos[index2-1]))
      tmpg = prevg + frac* (mapgpos[index2]-prevg)
      writegz([rs, pos, tmpg], writeout)
      index1 = index1+1
  elif pos > mappos[index2]:
    #current position in interpolation after marker
    if index2 == len(mappos)-1:
      #after the last site in the map (genetic position = maximum in map, note could try to extrapolate based on rate instead)
      writegz([rs, pos, mapgpos[index2]], writeout)
      index1 = index1+1
    else:
      #increment the marker
      index2 = index2+1

outfile.close()
