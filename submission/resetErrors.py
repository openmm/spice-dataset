from qcportal import PortalClient
from qcportal.record_models import RecordStatusEnum
import sys

dataset_name = sys.argv[1]
client = PortalClient.from_file()
ds = client.get_dataset('singlepoint', dataset_name)
ids = [r.id for e, s, r in ds.iterate_records(status=RecordStatusEnum.error)]
client.reset_records(ids)