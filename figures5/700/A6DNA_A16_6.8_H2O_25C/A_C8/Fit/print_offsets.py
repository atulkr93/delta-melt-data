import sys
import numpy as np
import pandas as pd

bf = pd.read_csv(str(sys.argv[1]))

slps = np.sort(np.unique(bf['SLP'].values))
onres_slps = []

for dummy_slp in slps:
    bf_subset = bf.loc[bf['SLP'] == dummy_slp]
    offsets = np.sort(np.array(bf_subset['Offset']*-1.0))
    offset_string = "{" + str(offsets[0])
    #print offsets 
    onres_flag = 0
    if len(offsets) == 1:
        onres_slps.append(dummy_slp)
        onres_flag = 1
    else:
        for dummy in range(1, len(offsets)):
            if offsets[dummy] != 0.0:
                offset_string = offset_string + ", " + str(offsets[dummy])
            else:
                onres_slps.append(dummy_slp)
    if onres_flag == 0:
        offset_string = offset_string + "}"
        print '[' + str(dummy_slp) + ']' + ", " + offset_string

print onres_slps, '{0}'
