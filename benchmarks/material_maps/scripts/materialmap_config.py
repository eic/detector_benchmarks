## For ACTS material map
## Read config-map.json and turn on mapping for approach 1 and/or 2 of each sensitive surface.
## Default binning is used
## Shujie Li, 10,2023

import pandas as pd
import numpy as np
import json
dirname = "./materialmap_no_beampipe/"
fname   = "config-map.json"
out_name = "config-map_new.json"

## turn on mapmaterial for all approach 1
## also turn on approach 2 for endcaps

f  = open(dirname+fname)
dd = json.load(f)

## get all endcap volume numbers
v_endcap = [] # list of volume numbers
n_endcap = [] # list of names
ee=dd['Volumes']['entries']
print ("Volume ID        Name       Approaches")
for vv in np.arange(len(ee)):
    nn = ee[vv]['value']['NAME']
    if  "|" not in nn and "Gap"  not in nn:
        if "Endcap" in nn:
            print(ee[vv]['volume'], nn, "1,2")#print(ee[vv]['value'])#['NAME'])
            v_endcap.append(ee[vv]['volume'])
            n_endcap.append(nn)
        else:
            print(ee[vv]['volume'], nn,"1")#print(ee[vv]['value'])#['NAME'])


for vv in np.arange(1,1+len(dd['Surfaces'])):
    for ii,tt in enumerate(dd['Surfaces'][str(vv)]):
        if 'approach' in tt:
            if tt['approach']==1 or (vv in v_endcap and tt['approach']==2):
                # print (vv,tt)
                # if "MPGD" not in tt:
                dd['Surfaces'][str(vv)][ii]['value']['material']['mapMaterial']=True
                # print(vv,dd['Surfaces'][str(vv)][ii])

with open(dirname+out_name, "w") as outfile:
    json.dump(dd, outfile, indent=4)
print("\n---------")
print("DONE! Output at ",dirname+out_name)
print("WARNING: you may still need to MANNUALLY turn on mapping and adjust binning for the cylindrical boundary of central beampipe, which is often vol 2, boundary 4, ")
print("")
