#!/usr/bin/env python3
"""
Collect runtimes from multiple log files and output a single row with user times.

Run like this:
runtimes.py $(ls -v phased/whatshap/realign/sim.pacbio.child.chr1.cov*.log)
"""
import sys
import re

def last_line(path):
	with open(path) as f:
		line = None
		for line in f:
			pass
	return line

common_prefix = sys.argv[1]


# example line:
# 208588 kB; real 440.70; user 347.76; sys 2.19
regex = re.compile('(?P<kb>[0-9]*) kB; real (?P<real>[0-9.]*); user (?P<user>[0-9.]*); sys (?P<sys>[0-9.]*)\n')

memory = []
user_time = []

for path in sys.argv[1:]:
	line = last_line(path)
	while not path.startswith(common_prefix):
		common_prefix = common_prefix[:-1]
	match = regex.match(line)
	memory.append(match.group('kb'))
	user_time.append(int(float(match.group('user'))+0.5))
print(common_prefix, *user_time, sep='\t')
