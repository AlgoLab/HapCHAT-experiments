#!/usr/bin/env python3
import re
from argparse import ArgumentParser
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def read_table(path):
	table = pd.read_csv(path, sep='\t')

	def nametocov(s):
		g = re.search('.cov([0-9]+|all)', s).group(1)
		if g == 'all':
			return 20
		else:
			return int(g)

	table['coverage'] = table['file_name1'].apply(nametocov)

	# maps dataset name to (algorithm, indels_or_not) pairs
	DS_ALGO = {
		'whatshap/classic': ('whatshap-norealign', False),
		'whatshap/trio': ('whatshap-trio', False),
		'whatshap/realign': ('whatshap', False),
		'read-backed-phasing': ('read-backed-phasing', False),
	}
	def algo_indels(s):
		if s in DS_ALGO:
			return DS_ALGO[s]
		fields = s.split('/')
		assert fields[1] in ('indels', 'noindels')
		return fields[0], fields[1] == 'indels'

	def simreal(s):
		if 'sim.child' in s:
			return 'sim'
		elif 'ashk.child' in s:
			return 'real'

	table['algorithm'] = table['dataset_name1'].apply(lambda ds: algo_indels(ds)[0])
	table['indels'] = table['dataset_name1'].apply(lambda ds: algo_indels(ds)[1])
	table['individual'] = table['file_name0'].apply(simreal)
	table.set_index(['indels', 'individual', 'algorithm', 'coverage'], inplace=True)
	table.sort_index(inplace=True)
	return table


def plot(fig, table, yfunc, yaxislabel, loc=0, yscale=None, ylim=None):
	ALGORITHMS = [
		('hapcut', 'hapCUT', 'o:'),
		('read-backed-phasing', 'ReadBackedPhasing', 'o-.'),
		('phaser', 'phASER', 'o--'),
		('whatshap-norealign', 'WhatsHap (no re-align)', 's:'),
		('whatshap', 'WhatsHap', 's-')
	]

	for individual, title, index in [('sim', 'Simulated data', 1), ('real', 'Real data', 2)]:
		ax = fig.add_subplot(1, 2, index)
		individual_table = table.xs(individual, level='individual')

		# Plot all curves
		for algorithm, name, marker in ALGORITHMS:
			t = individual_table.xs(algorithm, level='algorithm').reset_index()
			if t.empty:
				continue
			ax.plot(t['coverage'], yfunc(t), marker, label=name, markersize=5)
			print(yaxislabel, 'at coverage', list(t['coverage'])[-1], 'for tool', algorithm, 'is', list(yfunc(t))[-1])

		# Setup axes

		ax.set_title(title)
		if index == 1:  # y-axis label only on left plot
			ax.set_ylabel(yaxislabel)

		# To despine, we would need more horizontal space between subplots
		#sns.despine(ax=ax, right=False)
		ax.xaxis.set_ticks_position('bottom')
		ax.legend(loc=loc)

		ax.set_xticks([2, 3, 4, 5, 10, 15, 20])
		ax.set_xticklabels(['2', '3', '4', '5', '10', '15', 'all'])
		ax.set_xlabel('Coverage')
		ax.set_xlim(1, 21)

		if ylim:
			ax.set_ylim(ylim)
		if yscale:
			ax.set_yscale(yscale)
		ticklabels = ax.get_yticks()
		ax.set_yticklabels(['{:g}%'.format(y * 100) for y in ticklabels])


def main():
	parser = ArgumentParser()
	parser.add_argument('--type', choices=('switches', 'connections', 'indelswitches', 'indelphased'), default='switches')
	parser.add_argument('summary_eval', help='summary.eval file')
	parser.add_argument('plot', help='path to output PDF')
	args = parser.parse_args()

	table = read_table(args.summary_eval)
	indel_table = table.xs(True, level='indels').copy()
	noindel_table = table.xs(False, level='indels').copy()
	del table

	indel_table['extra_phased_snvs'] = (indel_table['phased_snvs'] - noindel_table['phased_snvs'])
	indel_table['extra_phased'] = (indel_table['phased'] - noindel_table['phased'])
	indel_table['extra_switches'] = indel_table['all_switches'] - noindel_table['all_switches']

	n_real_snvs = noindel_table.loc[('real', 'whatshap', 2)]['heterozygous_snvs']
	n_sim_snvs = noindel_table.loc[('sim', 'whatshap', 2)]['heterozygous_snvs']
	noindel_table.loc[('real',), 'connection_ratio'] = noindel_table.loc[:, 'all_assessed_pairs'] / (n_real_snvs - 1)
	noindel_table.loc[('sim',), 'connection_ratio'] = noindel_table.loc[:, 'all_assessed_pairs'] / (n_sim_snvs - 1)

	sns.set_style('ticks')
	pal = sns.color_palette()
	pal = pal[2:] + pal[:2]
	with sns.color_palette(pal):
		fig = plt.figure(figsize=(8.5, 5))

		if args.type == 'switches':
			plot(fig, noindel_table, lambda t: t['all_switch_rate'], yaxislabel='Switch error rate', yscale='log', ylim=(0.0001, 0.18))
		elif args.type == 'connections':
			plot(fig, noindel_table, lambda t: t['connection_ratio'], yaxislabel='Phase connection ratio', loc='lower right')
		elif args.type == 'indelswitches':
			plot(fig, indel_table, lambda t: t['extra_switches'] / t['extra_phased'], 'Extra switches per extra phased variants')
		elif args.type == 'indelphased':
			func = lambda t: (t['phased'] - t['phased_snvs']) / (t['heterozygous_variants'] - t['heterozygous_snvs'])
			plot(fig, indel_table, func, yaxislabel='Phased non-SNVs')

		fig.tight_layout()
		fig.savefig(args.plot)


if __name__ == '__main__':
	main()

