from qcportal import PortalClient
from collections import defaultdict
import sys

dataset_name = sys.argv[1]
client = PortalClient('https://ml.qcarchive.molssi.org')
ds = client.get_dataset('singlepoint', dataset_name)
status = defaultdict(int)
for e, s, r in ds.iterate_records():
    status[r.status] += 1
for key in status:
    print(key, status[key])