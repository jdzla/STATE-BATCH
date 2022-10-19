import os, sys, json, subprocess
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-i', '--run_idx', type=int, default=1)
args = parser.parse_args()
index = args.run_idx

with open('atom_run_jobdict.json', 'r') as read_file:
  jobdict = json.load(read_file)

os.chdir(jobdict[index]['dirname'])
subprocess.run(jobdict[index]['command'], shell=True)
