from qcportal import PortalClient
import sys

dataset_name = sys.argv[1]
client = PortalClient('https://ml.qcarchive.molssi.org')
ds = client.get_dataset('singlepoint', dataset_name)
ds.print_status()
