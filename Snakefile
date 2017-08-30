# kate: syntax Python; space-indent off; indent-mode python; indent-width 4;
"""
How to re-run this workflow
---------------------------


Install dependencies:
- (see Dockerfile)
- make 'gatk' binary available
- make 'shapeit' binary available
- copy human reference and BWA index into reference/ if you already have it
  (downloaded and generated otherwise)

sudo docker build -t whatshap-experiments docker/
sudo docker run -it -v $PWD:/io/ whatshap-experiments snakemake -np
"""

import pysam
import textwrap

shell.executable("/bin/bash")

picard_tmpdir_switch = ''
if 'TMPDIR' in os.environ:
	picard_tmpdir_switch = 'TMP_DIR=%s' % os.environ['TMPDIR']

datasets = ['pacbio', 'ont']
coverage = [2, 3, 4, 5, 10, 15, 'all']
reference = 'reference/GRCh38_full_analysis_set_plus_decoy_hla.fa'

time = "/usr/bin/time -f '%M kB; real %e; user %U; sys %S'"

# Paths to scripts distributed along with the Snakefile
eval_plot = srcdir('scripts/evalplot.py')

# Tools assumed to be installed somewhere on the PATH.
phaser = 'phaser'
extract_hairs = 'extractHAIRS'
hapcut = 'HAPCUT'
extract_hairs2 = 'extractHAIRS2'  # this is the extractHAIRS that comes with hapCUT 2
hapcut2 = 'HAPCUT2'
whatshap = 'whatshap'

PROGRAMS = [
	'hapcut/indels',
	'hapcut/noindels',
	'hapcut2/indels',
	'hapcut2/noindels',
	'read-backed-phasing/noindels',
	'phaser/indels',
	'phaser/noindels',
	'whatshap-norealign/noindels',
	'whatshap/noindels',
	'whatshap/indels',
]

# Software that must be installed manually prior to running
# the Snakefile due to licensing restrictions
gatk_jar = 'restricted-software/GenomeAnalysisTK.jar'

dataset_pattern = '{dataset,[a-z]+}.chr{chromosome,[0-9]+}'


rule master:
	input:
		'reference/OK',
		'eval/summary.eval',
		expand('eval/{type}.pdf', type=['switches', 'connections', 'indelswitches', 'indelphased'])


rule clean:
	shell:
		# download/ not included
		"rm -rf bam/ stats/ vcf/ phased/ eval/ hairs/"


## Rules for external dependencies: Software and data files

rule gatk_missing:
	output: gatk_jar
	shell:
		"""
		echo
		echo "Due to licensing restrictions, you need to manually download"
		echo "GATK 3.5 (3.5-0-g36282e4) and make it available to this pipeline."
		echo
		echo "Choose one of the following options:"
		echo
		echo " - If you run from the Docker image, then download GATK 3.5,"
		echo "   unpack it, place the jar file at '{output}'"
		echo "   and then re-run snakemake."
		echo
		exit 1
		"""


rule references_ok:
	output:
		'reference/OK'
	input:
		expand('reference/GRCh38_full_analysis_set_plus_decoy_hla.{ext}',
			ext=['dict', 'fa', 'fa.alt', 'fa.amb', 'fa.ann', 'fa.bwt', 'fa.fai', 'fa.pac', 'fa.sa'])
	shell:
		"cd reference && md5sum -c MD5SUM && touch OK"


rule download_reference:
	output:
		'reference/GRCh38_full_analysis_set_plus_decoy_hla.{ext}'
	shell:
		"""
		wget -O {output}.incomplete ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.{wildcards.ext}
		mv {output}.incomplete {output}
		"""


rule download_platinum:
	output:
		'platinum/NA12878.vcf.gz'
	shell:
		"""
		wget --ftp-password= --ftp-user=platgene_ro -O {output}.incomplete ftp://ussd-ftp.illumina.com/2017-1.0/hg38/small_variants/NA12878/NA12878.vcf.gz
		mv {output}.incomplete {output}
		cd platinum && md5sum -c MD5SUM
		"""


rule download_pacbio:
	output:
		'pacbio/hg38.NA12878-WashU.chr1.bam'
	threads: 4
	shell:
		# The conda samtools packages do not support HTTPS URLs, so use the
		# system samtools if it is available
		"""
		if test -x /usr/bin/samtools; then SAMTOOLS=/usr/bin/samtools; else SAMTOOLS=samtools; fi
		$SAMTOOLS view -@ 4 -b https://downloads.pacbcloud.com/public/dataset/na12878/hg38.NA12878-WashU.bam chr1 > {output}
		# TODO will MD5SUM always be the same?
		"""


rule download_nanopore:
	output:
		'nanopore/chr1.sorted.bam'
	shell:
		"""
		wget -O {output}.incomplete http://s3.amazonaws.com/nanopore-human-wgs/chr1.sorted.bam
		mv {output}.incomplete {output}
		cd nanopore && md5sum -c MD5SUM
		"""


rule symlink_bam:
	output:
		pacbio='bam/pacbio.chr1.covall.bam',
		ont='bam/ont.chr1.covall.bam'
	input:
		pacbio='pacbio/hg38.NA12878-WashU.chr1.bam',
		ont='nanopore/chr1.sorted.bam'
	shell:
		"""
		ln -srf {input.pacbio} {output.pacbio}
		ln -srf {input.ont} {output.ont}
		"""


rule calculate_coverage:
	input: 'bam/' + dataset_pattern + '.cov{coverage,(all|[0-9]+)}.bam'
	output: 'stats/bam/' + dataset_pattern + '.cov{coverage,(all|[0-9]+)}.coverage'
	message: 'Computing coverage for {input}'
	run: 
		bam = pysam.Samfile(input[0])
		length = None
		for e in bam.header.get('SQ'):
			if e['SN'] == wildcards.chromosome:
				length = e['LN']
		assert length != None
		shell("samtools depth {input} | awk '{{sum+=$3}} END {{ print sum/{length} }}' > {output}")


rule downsample:
	input:
		bam='bam/' + dataset_pattern + '.covall.bam',
		coverage='stats/bam/' + dataset_pattern + '.covall.coverage'
	output: 
		bam='bam/' + dataset_pattern + '.cov{coverage,([0-9]+)}.bam',
		bai='bam/' + dataset_pattern + '.cov{coverage,([0-9]+)}.bai'
	log: 'bam/' + dataset_pattern + '.cov{coverage,([0-9]+)}.log'
	message: 'Downsampling {input.bam} to {wildcards.coverage}x'
	run:
		input_coverage = float(open(input.coverage).readline())
		p = float(wildcards.coverage) / input_coverage
		seed = hash(output)
		shell('picard DownsampleSam INPUT={input.bam} RANDOM_SEED={seed} CREATE_INDEX=true OUTPUT={output.bam} PROBABILITY={p} VALIDATION_STRINGENCY=SILENT > {log} 2>&1')


rule unphase:
	input: '{base}.phased.vcf'
	output: '{base}.unphased.vcf'
	shell: '{whatshap} unphase {input} > {output}'


## Rules for running the phasing tools.
## Some have more than one rule if pre- and/or postprocessing is needed

rule gatk_read_backed_phasing:
	input:
		bam='bam/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.bam',
		bai='bam/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.bai',
		vcf='vcf/unphased.chr{chromosome,[0-9]+}.vcf',
		ref=reference,
		dictfile=reference.replace('.fasta', '.dict'),
		gatk_jar=gatk_jar
	output:
		vcf='phased/read-backed-phasing/noindels/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.vcf',
	log: 'phased/read-backed-phasing/noindels/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.log'
	message: 'Running GATK\'s ReadBackedPhasing on {input.bam}'
	threads: 1  # ReadBackedPhasing cannot run with more than one thread
	run:
		#input_line = ' '.join('-I ' + bam for bam in input.bams)
		shell(r"""
		{time} java -jar {gatk_jar} -T ReadBackedPhasing \
			-R {input.ref} \
			-I {input.bam} \
			-L {input.vcf} \
			-V {input.vcf} \
			--phaseQualityThresh 1 --min_base_quality_score 1 \
			-o {output.vcf} >& {log}""")


hapcut_out = dataset_pattern + '.cov{coverage,([0-9]+|all)}'


rule extract_hairs:
	"""hapCUT-specific pre-processing"""
	output:
		txt='hairs/{indelsornot,(indels|noindels)}/' + hapcut_out + '.txt'
	input:
		ref=reference,
		bam='bam/{dataset}.chr{chromosome}.cov{coverage}.bam',
		bai='bam/{dataset}.chr{chromosome}.cov{coverage}.bai',
		vcf='vcf/unphased.chr{chromosome,[0-9]+}.vcf'
	log:
		'hairs/{indelsornot,(indels|noindels)}/' + hapcut_out + '.log'
	run:
		extra = ' --indels 1' if wildcards.indelsornot == 'indels' else ''
		shell("{time} {extract_hairs} --ref {input.ref} {extra} --VCF {input.vcf} --bam {input.bam} --maxIS 600 > {output.txt} 2> {log}")


rule hapcut:
	output:
		txt='phased/hapcut/{indelsornot,(indels|noindels)}/' + hapcut_out + '.txt'
	input:
		txt='hairs/{indelsornot}/' + hapcut_out + '.txt',
		vcf='vcf/unphased.chr{chromosome,[0-9]+}.vcf'
	log: 'phased/hapcut/{indelsornot}/' + hapcut_out + '.phase.log'
	shell:
		"{time} {hapcut} --fragments {input.txt} --VCF {input.vcf} --output {output.txt} >& {log}"


# We use the extractHAIRS distributed with hapCUT as the one distributed with
# hapCUT 2 gives unusable output (perhaps the same problem as the one that
# we fixed in the hapCUT 1 version with SAM = and X operators)
rule hapcut2:
	output:
		txt='phased/hapcut2/{indelsornot,(indels|noindels)}/' + hapcut_out + '.txt'
	input:
		txt='hairs/{indelsornot}/' + hapcut_out + '.txt',
		vcf='vcf/unphased.chr{chromosome,[0-9]+}.vcf',
	log: 'phased/hapcut2/{indelsornot}/' + hapcut_out + '.phase.log'
	shell:
		"{time} {hapcut2} --fragments {input.txt} --VCF {input.vcf} --output {output.txt} >& {log}"


rule hapcut_to_vcf:
	output:
		vcf='phased/hapcut/{indelsornot,(indels|noindels)}/' + hapcut_out + '.vcf'
	input:
		vcf='vcf/unphased.chr{chromosome,[0-9]+}.vcf',
		hapcut='phased/hapcut/{indelsornot}/' + hapcut_out + '.txt'
	log: 'phased/hapcut/{indelsornot}/' + hapcut_out + '.vcf.log'
	shell:
		"{time} {whatshap} hapcut2vcf {input.vcf} {input.hapcut} > {output.vcf} 2> {log}"


rule hapcut2_to_vcf:
	output:
		vcf='phased/hapcut2/{indelsornot,(indels|noindels)}/' + hapcut_out + '.vcf'
	input:
		vcf='vcf/unphased.chr{chromosome,[0-9]+}.vcf',
		hapcut='phased/hapcut2/{indelsornot}/' + hapcut_out + '.txt'
	log: 'phased/hapcut2/{indelsornot}/' + hapcut_out + '.vcf.log'
	shell:
		"{time} {whatshap} hapcut2vcf {input.vcf} {input.hapcut} > {output.vcf} 2> {log}"


rule phaser:
	output:
		vcf='phased/phaser/{indelsornot,(indels|noindels)}/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.vcf'
	params:
		base='phased/phaser/{indelsornot}/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}'
	input:
		bam='bam/{dataset}.chr{chromosome}.cov{coverage}.bam',
		bai='bam/{dataset}.chr{chromosome}.cov{coverage}.bai',
		vcf='vcf/unphased.chr{chromosome,[0-9]+}.vcf.gz',
		tbi='vcf/unphased.chr{chromosome,[0-9]+}.vcf.gz.tbi',
	log:
		'phased/phaser/{indelsornot}/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.log'
	run:
		extra = ' --include_indels 1' if wildcards.indelsornot == 'indels' else ' --include_indels 0'
		sample = 'NA12878'
		shell("{time} {phaser}{extra} --bam {input.bam} --write_vcf 1 --gw_phase_vcf 2 --pass_only 0 --vcf {input.vcf} --sample {sample} --mapq 1 --baseq 1 --paired_end 0 --o {params.base} >& {log}")
		shell("gunzip -f {params.base}.vcf.gz")
		#shell(r"sed -i.orig '/^#/!s|\bGT:PG:|PG:GT:|;/^#/!s|:PI:|:PS:|' {params.base}.vcf")


rule whatshap_norealign:
	input:
		bam='bam/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.bam',
		vcf='vcf/unphased.chr{chromosome,[0-9]+}.vcf',
	output: 'phased/whatshap-norealign/noindels/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.vcf'
	log: 'phased/whatshap-norealign/noindels/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.log'
	shell: '{time} {whatshap} phase {input.vcf} {input.bam} > {output} 2> {log}'


rule whatshap_noindels:  # with re-alignment
	input:
		ref=reference,
		bam='bam/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.bam',
		vcf='vcf/unphased.chr{chromosome,[0-9]+}.vcf',
	output:
		vcf='phased/whatshap/noindels/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.vcf',
		log='phased/whatshap/noindels/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.log'
	shell: '{time} {whatshap} phase --reference {input.ref} {input.vcf} {input.bam} > {output.vcf} 2> {output.log}'


rule whatshap_indels:  # with re-alignment
	input:
		ref=reference,
		bam='bam/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.bam',
		vcf='vcf/unphased.chr{chromosome,[0-9]+}.vcf',
	output:
		vcf='phased/whatshap/indels/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.vcf',
		log='phased/whatshap/indels/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.log'
	shell: '{time} {whatshap} phase --indels --reference {input.ref} {input.vcf} {input.bam} > {output.vcf} 2> {output.log}'


## Evaluation: compare phasing results to ground truth, get VCF statistics


rule evaluate_phasing_tool:
	input:
		truth='vcf/truth.chr{chromosome}.vcf',
		phased='phased/{program}/' + dataset_pattern + '.cov{coverage}.vcf'
	output:
		'eval/{program,(whatshap-norealign/noindels|whatshap/trio|whatshap/noindels|whatshap/indels|read-backed-phasing/noindels|hapcut/indels|hapcut/noindels|hapcut2/indels|hapcut2/noindels|phaser/indels|phaser/noindels)}/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.eval'
	log:
		'eval/{program}/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.log'
	shell:
		'{whatshap} compare --names truth,{wildcards.program} --tsv-pairwise {output} {input.truth} {input.phased} 2>&1 > {log}'


rule whatshap_stats:
	output:
		stats='eval/{program}/' + dataset_pattern + '.cov{coverage}.stats'
	input:
		vcf='phased/{program}/' + dataset_pattern + '.cov{coverage}.vcf'
	log:
		'eval/{program}/' + dataset_pattern + '.cov{coverage}.statslog'
	shell:
		'{whatshap} stats --tsv {output.stats} {input.vcf} &> {log}'


rule connected_components:
	# Extract the number of connected components from WhatsHap log output
	output:
		components='eval/components/{indelsornot,(indels|noindels)}/' + dataset_pattern + '.cov{coverage}.components'
	input:
		whatshap_log='phased/whatshap/{indelsornot}/' + dataset_pattern + '.cov{coverage}.log'
	run:
		r = re.compile(r'Best-case phasing would result in [0-9]+ non-singleton phased blocks \((?P<blocks>[0-9]+) in total\).*')
		with open(input.whatshap_log) as f:
			components = None
			for line in f:
				m = r.match(line)
				if m:
					components = int(m.group('blocks'))
					# no break - we want the last occurrence
			if components is None:
				sys.exit("WhatsHap log file missing line with 'Best-case phasing would result ...'")
		with open(output.components, 'w') as f:
			print('dataset', 'indels', 'coverage', 'components', sep='\t', file=f)
			indels = int(wildcards.indelsornot == 'indels')
			print(wildcards.dataset, indels, wildcards.coverage, components, sep='\t', file=f)



## Create two tables with evaluation results

def all_eval_paths(extension):
	l1 = expand('eval/{program}/{dataset}.chr1.cov{coverage}' + extension,
		program=PROGRAMS,
		dataset=datasets, coverage=coverage),
	# Exclude because they did not finish within 24h
	#l2 = tuple(s for s in l2 if not ('/ashk.' in s and '.covall.' in s))
	return l1 + l2


rule evaluation_summary:
	input:
		evals=all_eval_paths('.eval'),
		stats=all_eval_paths('.stats')
	output: 'eval/summary.eval'
	shell:
		"""
		(
			paste {input.evals[0]} {input.stats[0]} | head -n1
			for e in {input.evals}; do
				s=${{e%%.eval}}.stats
				paste $e $s | sed 1d
			done
		) > {output}
		"""


rule summarize_connected_components:
	output:
		summary='eval/summary.components'
	input:
		components=expand('eval/components/{indels}/{dataset}.chr1.cov{coverage}.components',
			indels=['indels', 'noindels'], dataset=datasets, coverage=coverage)
	shell:
		"""
		(
			head -n1 {input.components[0]}
			for f in {input.components}; do
				sed 1d $f
			done
		) > {output}
		"""


rule plot_eval:
	output:
		pdf='eval/{type}.pdf'
	input:
		eval='eval/summary.eval',
		components='eval/summary.components'
	run:
		shell("{eval_plot} --type={wildcards.type} {input.eval} {input.components} {output.pdf}")


## General rules

rule vcf_tabix:
	input: '{path}.vcf.gz'
	output: '{path}.vcf.gz.tbi'
	shell: 'tabix {input}'


rule vcf_bgzip:
	input: '{path}.vcf'
	output: '{path}.vcf.gz'
	shell: 'bgzip < {input} > {output}'


#rule index_reference:
	#output:
		#'{path}.fasta.bwt'
	#input:
		#'{path}.fasta'
	#shell:
		#"bwa index {input}"


#rule CreateSequenceDict:
	#output: '{base}.dict'
	#input: '{base}.fasta'
	#shell:
		#"picard CreateSequenceDictionary R={input} O={output}"
