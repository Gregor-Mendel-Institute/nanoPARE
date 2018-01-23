import argparse

########################
### ARGUMENT PARSING ###
########################

parser = argparse.ArgumentParser()
parser.add_argument('-R','--rpm',dest='RPM',
                    help="Path to a RPM table.",
                    required=True)
parser.add_argument('-T','--tpm',dest='TPM',
                    help="Path to TPM table.",
                    required=True)
parser.add_argument('-W','--winsorize',dest='WINSORIZE',
                    help="Exclude outliers by examining the middle 90% of the data.",
                    default=False, action='store_true')
args = parser.parse_args()

rpm_vals = {}
tpm_vals = {}

for line in open(args.RPM):
    l = line.rstrip().split()
    rpm_vals[l[0]] = float(l[1])

for line in open(args.TPM):
    l = line.rstrip().split()
    tpm_vals[l[0]] = float(l[1])

ratios = []
rpm_ids = set(list(rpm_vals.keys()))
tpm_ids = set(list(tpm_vals.keys()))

shared_ids = list(rpm_ids.intersection(tpm_ids))
for i in shared_ids:
    if rpm_vals[i] > 0 and tpm_vals[i] > 0:
        ratios += [rpm_vals[i]/tpm_vals[i]]

if args.WINSORIZE:
    ratios.sort()
    threshold = int(len(ratios)*0.05)
    ratios = ratios[threshold:-threshold]

print(float(sum(ratios))/len(ratios))
