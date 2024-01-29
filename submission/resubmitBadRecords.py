from qcportal import PortalClient
import numpy as np
import sys

# Look for suspect records that should be rerun.

dataset_name = sys.argv[1]
client = PortalClient.from_file()
dataset = client.get_dataset('singlepoint', dataset_name)
include_specs = [name for name in dataset.specifications if dataset.specifications[name].specification.method == 'wb97m-d3bj']
recs = list(dataset.iterate_records(specification_names=include_specs))
specifications = list(set(s for e, s, r in recs))
bad_recs = []
for e, s, r in recs:
    if r is not None and r.status == 'complete':
        if np.max(np.abs(r.properties['dft total gradient'])) > 1 or r.compute_history[-1].provenance.version == '1.4a3.dev63':
            bad_recs.append((e, s, r))
if len(bad_recs) == 0:
    print('All records look ok')
    exit()
print('Found', len(bad_recs), 'suspect records')

# Create the new dataset by cloning the old one.

new_name = f'{dataset_name[:-1]}{int(dataset_name[-1])+1}'
print('Creating new dataset:', new_name)
new_dataset = client.add_dataset('singlepoint', new_name)
print('ID:', new_dataset.id)
new_dataset.add_entries(list(dataset.iterate_entries()))
for name in specifications:
    new_dataset.add_specification(name, dataset.specifications[name].specification)
new_dataset.submit()

# Verify that all records were matched correctly.

old_status = dataset.status()
new_status = new_dataset.status()
for spec in specifications:
    if old_status[spec] != new_status[spec]:
        print('Status is not correct.  Original:\n')
        dataset.print_status()
        print('\nNew:\n')
        new_dataset.print_status()
        exit()

# Remove the records we want to rerun and then submit the dataset again.

rerun_entries = list(set(e for e, s, r in bad_recs))
rerun_specifications = list(set(s for e, s, r in bad_recs))
new_dataset.invalidate_records(rerun_entries, rerun_specifications)
new_dataset.remove_records(rerun_entries, rerun_specifications, delete_records=False)
new_dataset.submit(rerun_entries, rerun_specifications, find_existing=False, tag='spice-psi4-181-2')
print('Submission complete\n')
new_dataset.print_status()
