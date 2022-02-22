import numpy as np
import sys
import pandas as pd
bf = pd.read_csv(sys.argv[1])

for dummy in np.unique(bf['SLP'].values):
    bf_new = bf[bf['SLP'] == dummy]
    print dummy, np.min(bf_new['Offset'].values), np.max(bf_new['Offset'].values)

