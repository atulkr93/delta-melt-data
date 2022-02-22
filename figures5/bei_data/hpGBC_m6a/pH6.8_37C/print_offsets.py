import sys
import numpy as np
import pandas as pd

bf = pd.read_csv(str(sys.argv[1]))

slps = np.sort(np.unique(bf['slp(hz)'].values))

for dummy_slp in slps:
    bf_subset = bf.loc[bf['slp(hz)'] == dummy_slp]
    offsets = np.sort(np.array(bf_subset['offset(hz)']))
    offset_string = "{" + str(offsets[0]) 
    for dummy in range(1, len(offsets)):
        offset_string = offset_string + ", " + '%.1f'%offsets[dummy]
    offset_string = offset_string + "}"
    print '[' + '%.1f'%dummy_slp + ']' + ", " + offset_string

