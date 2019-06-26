from RAW_create import crete_raw
from FC_create import create_fc

import time
start_time = time.time()

filename = "Moon"
voxel_size = 0.1
layers = -1

raw_filename = crete_raw(filename=filename, voxel_size=voxel_size, layers=layers)
create_fc(raw_filename)
print("--- %s seconds ---" % (time.time() - start_time))

