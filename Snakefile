# kate: syntax Python; space-indent off; indent-mode python; indent-width 4;
"""
data
ONT: 
ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/
https://github.com/nanopore-wgs-consortium/NA12878

PacBio : https://downloads.pacbcloud.com/public/dataset/na12878/
 
How to re-run this workflow
---------------------------


Install dependencies:
- (see Dockerfile)
- make 'gatk' binary available
- make 'shapeit' binary available
- copy human reference and BWA index into reference/ if you already have it
  (downloaded and generated otherwise)

cd docker
sudo docker build -t whatshap-experiments .

docker run -it -v $PWD:/io/ whatshap-experiments snakemake -np


"""

import pysam
import textwrap

picard_tmpdir_switch = ''
if 'TMPDIR' in os.environ:
	picard_tmpdir_switch = 'TMP_DIR=%s' % os.environ['TMPDIR']

datasets = ['ceph']
individuals = ['child']  #['mother', 'father', 'child']
coverage = [2, 3, 4, 5, 10, 15, 'all']
platform = ['oxford', 'pacbio']
platform_to_ref = {'pacbio': 'reference/hg38.chr1.fa', 'oxford': 'reference/GRCh38_full_analysis_set_plus_decoy_hla.chr1.fa'}
reference = 'reference/human_g1k_v37.fasta'
role_to_sampleid = {'child':'NA12878'}

time = "/usr/bin/time -f '%M kB; real %e; user %U; sys %S'"

# Paths to scripts distributed along with the Snakefile
eval_plot = srcdir('scripts/evalplot.py')
vcf_merge_trio = srcdir('scripts/vcf_merge_trio.pl')
genomesimulator = srcdir('scripts/genomesimulator.py')
artificial_child = srcdir('scripts/artificial-child.py')

# Tools assumed to be installed somewhere on the PATH.
phaser = 'phaser'
extract_hairs = 'extractHAIRS'
hapcut = 'HAPCUT'
extract_hairs2 = 'extractHAIRS2'  # this is the extractHAIRS that comes with hapCUT 2
hapcut2 = 'HAPCUT2'
gzip = 'pigz'
whatshap = 'whatshap'

# Software that must be installed manually prior to running
# the Snakefile due to licensing restrictions
shapeit = 'restricted-software/shapeit'
gatk_jar = 'restricted-software/GenomeAnalysisTK.jar'

# hapcompass is unused
hapcompass = 'java -Xmx220g -jar .../hapcompass-0.8.1.jar'

dataset_pattern = '{dataset,[a-z]+}.{platform,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}'


rule master:
	input:
		'eval/summary.eval',
		expand('eval/{type}.pdf', type=['switches', 'connections', 'indelswitches', 'indelphased'])


rule clean:
	shell:
		# download/ not included
		"rm -rf bam/ fastq/ genmap/ shapeit/ sim/ stats/ vcf/ phased/ eval/ hairs/"


## Rules for external dependencies: Software and data files

rule shapeit_missing:
	output: shapeit
	shell:
		"""
		echo
		echo "Due to licensing restrictions, you need to manually download "
		echo "shapeit and make it available to this pipeline."
		echo
		echo "Go to https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#download ,"
		echo "download version v2.r837 (shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz),"
		echo "unpack the tgz file and place the bin/shapeit binary into '{shapeit}'."
		echo "Then re-run snakemake."
		echo
		exit 1
		"""


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

#TODO: update link to point to new dataset
rule download_ashkenazim:
	threads: 100
	output:
		protected("download/AshkenazimTrio/{file}.{ext,(bam|bam.bai|vcf.gz)}")
	shell:
		"wget -O {output} ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/{wildcards.file}.{wildcards.ext}"


#TODO: update link to point to new dataset
rule download_reference:
	output:
		'reference/human_g1k_v37.fasta'
	shell:
		"""
		wget -O {output}.gz.incomplete ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
		mv {output}.gz.incomplete {output}.gz
		# gunzip fails with "decompression OK, trailing garbage ignored" because
		# the file is razf-compressed (gzip-compatible, but with an index at end)
		gunzip -f {output}.gz || true
		cd reference && md5sum -c MD5SUM
		"""


## Downsample BAM files, add read group information, etc.

#TODO: update link to point to new dataset
def ashkenazim_trio_bam(wildcards):
	individual = {
		'mother': dict(name='mother', id='4', na='24143'),
		'father': dict(name='father', id='3', na='24149'),
		'child': dict(name='son', id='2', na='24385'),
	}[wildcards.individual]
	individual['chromosome'] = wildcards.chromosome
	individual['ext'] = wildcards.ext
	return 'download/AshkenazimTrio/HG00{id}_NA{na}_{name}/PacBio_MtSinai_NIST/MtSinai_'\
			'blasr_bam_GRCh37/hg00{id}_gr37_{chromosome}.{ext}'.format(**individual)

#TODO: update link to point to new dataset
rule create_ashk_pacbio_links:
	input: ashkenazim_trio_bam
	output: 'bam/incorrect-readgroup/ashk.pacbio.{individual}.chr{chromosome,[0-9]+}.covall.{ext,(bam|bam.bai)}'
	shell: 'ln -fsrv {input} {output}'


rule calculate_coverage:
	input: 'bam/' + dataset_pattern + '.cov{coverage,(all|[0-9]+)}.bam'
	output: 'stats/bam/' + dataset_pattern + '.cov{coverage,(all|[0-9]+)}.coverage'
	message: 'Computing coverage for {input}'
	run: 
		bam = pysam.Samfile(input[0])
		length = None
		for e in bam.header.get('SQ'):
			if e['SN'] == "chr"+ wildcards.chromosome:
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


rule add_read_groups:
	"""Add RG, prepend sample names to read names, strip unused tags"""
	input: 'bam/incorrect-readgroup/ashk.{platform,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.cov{coverage,([0-9]+|all)}.bam'
	output: 
		bam='bam/ashk.{platform,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.cov{coverage,([0-9]+|all)}.bam',
		bai='bam/ashk.{platform,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.cov{coverage,([0-9]+|all)}.bai'
	log: 'bam/ashk.{platform,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.cov{coverage,([0-9]+|all)}.readgroup.log'
	run: 
		sample = role_to_sampleid[wildcards.individual]
		shell(textwrap.dedent(
			"""
			( samtools view -h {input} | \\
				sed '/^[^@]/ s|^|{wildcards.individual}_|g' | \\
				sed -re 's_\\t(dq|dt|ip|iq|st|sq|mq):Z:[^\\t]*__g' | \\
				samtools view -u - | \\
				picard AddOrReplaceReadGroups {picard_tmpdir_switch} CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT I=/dev/stdin O={output.bam} ID={wildcards.individual} LB=UNKNOWN PL=PACBIO PU=UNKNOWN SM={sample}) > {log} 2>&1
			"""))


rule bam_to_fastq:
	input: 'bam/{file}.bam'
	output: 'fastq/{file}.fastq'
	shell: 'bedtools bamtofastq -i {input} -fq {output}'


## Prepare the downloaded VCF files: Extract chromosome 1, filter, etc.


rule unzip_vcf:
	input: 'download/AshkenazimTrio/analysis/NIST_CallsIn2Technologies_05182015/HG{id}-multiall-fullcombine.vcf.gz'
	output: 'unphased-ashk/HG{id,[0-9]+}-multiall-fullcombine.vcf'
	shell: 'zcat {input} > {output}'


rule filter_vcfs:
	input:
		mother='unphased-ashk/HG004-multiall-fullcombine.vcf',
		father='unphased-ashk/HG003-multiall-fullcombine.vcf',
		child='unphased-ashk/HG002-multiall-fullcombine.vcf'
	output: 'vcf/ashk.trio.unphased.vcf'
	log: 'vcf/ashk.trio.unphased.vcf.log'
	message: 'Filtering and merging input VCFs to create {output}'
	version: 2
	shell:
		r"{vcf_merge_trio} {input.mother} {input.father} {input.child} | sed '/^[^#]/ s| |\t|g; s_0|1_0/1_g; s_1|0_0/1_g; s_1/0_0/1_g; s_1|1_1/1_g; s_0|0_0/0_g' > {output} 2> {log}"


rule split_vcf:
	input: 'vcf/NA12878.vcf'
	output: 'vcf/ceph.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.phased.vcf'
	message: 'Extracting chromosome {wildcards.chromosome} from {input}'
	shell: """awk '/^#/ || ($1 == "chr{wildcards.chromosome}")' {input} > {output}"""


#rule shapeit_create_ped:
	#input: 'vcf/ashk.trio.chr{chromosome,[0-9]+}.unphased.vcf'
	#output: 'shapeit/trio.chr{chromosome,[0-9]+}.ped'
	#log: 'shapeit/trio.chr{chromosome,[0-9]+}.ped.log'
	#message: 'Creating PED file {output}'
	#shell: 'vcf-to-ped.py {input} {output} > {log} 2>&1'


rule unphase:
	input: '{base}.phased.vcf'
	output: '{base}.unphased.vcf'
	shell: '{whatshap} unphase {input} > {output}'


## Compute "ground truth" phasing with SHAPEIT

rule download_1000GP:
	output: 'download/1000GP_Phase3/{filename}'
	shell: 'wget -O {output} http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/{wildcards.filename}'


rule shapeit_check:
	input:
		legend='download/1000GP_Phase3/1000GP_Phase3_chr{chromosome,[0-9]+}.legend.gz',
		genmap='download/1000GP_Phase3/genetic_map_chr{chromosome,[0-9]+}_combined_b37.txt',
		refhaps='download/1000GP_Phase3/1000GP_Phase3_chr{chromosome,[0-9]+}.hap.gz',
		refsamples='download/1000GP_Phase3/1000GP_Phase3.sample',
		vcf='vcf/ashk.trio.chr{chromosome,[0-9]+}.unphased.vcf',
		shapeit=shapeit,
	output: 'shapeit/trio.chr{chromosome,[0-9]+}.snp.strand', 'shapeit/trio.chr{chromosome,[0-9]+}.snp.strand.exclude'
	log: 'shapeit/trio.chr{chromosome,[0-9]+}.check.log'
	shell: '({shapeit} -check -V {input.vcf} -M {input.genmap} --input-ref {input.refhaps} {input.legend} {input.refsamples} --output-log shapeit/trio.chr{wildcards.chromosome} || true) > {log} 2>&1'


rule shapeit:
	input:
		legend='download/1000GP_Phase3/1000GP_Phase3_chr{chromosome,[0-9]+}.legend.gz',
		genmap='download/1000GP_Phase3/genetic_map_chr{chromosome,[0-9]+}_combined_b37.txt',
		refhaps='download/1000GP_Phase3/1000GP_Phase3_chr{chromosome,[0-9]+}.hap.gz',
		refsamples='download/1000GP_Phase3/1000GP_Phase3.sample',
		vcf='vcf/ashk.trio.chr{chromosome,[0-9]+}.unphased.vcf',
		exclude='shapeit/trio.chr{chromosome,[0-9]+}.snp.strand.exclude',
		shapeit=shapeit,
	output: 'shapeit/trio.chr{chromosome,[0-9]+}.haps', 'shapeit/trio.chr{chromosome,[0-9]+}.sample' 
	log: 'shapeit/trio.chr{chromosome,[0-9]+}.run.log'
	message: 'Running SHAPEIT on chromosome {wildcards.chromosome}'
	shell: '{shapeit} -V {input.vcf} --exclude-snp {input.exclude} -M {input.genmap} --input-ref {input.refhaps} {input.legend} {input.refsamples} -O shapeit/trio.chr{wildcards.chromosome} > {log} 2>&1'

# We don't need rules about shapeit and trio because we decided to consider only real data for single individual.
rule shapeit_convert_to_vcf:
	input:
		haps='shapeit/trio.chr{chromosome,[0-9]+}.haps',
		shapeit=shapeit,
	output:
		vcf='shapeit/trio.chr{chromosome,[0-9]+}.phased.vcf'
	log: 'shapeit/trio.chr{chromosome,[0-9]+}.phased.vcf.log'
	message: 'Converting SHAPEIT output for chromosome {wildcards.chromosome} to VCF'
	shell: '{shapeit} -convert --input-haps shapeit/trio.chr{wildcards.chromosome} --output-vcf {output.vcf} > {log} 2>&1'


rule split_shapeit_vcf:
	input: 'shapeit/trio.chr{chromosome,[0-9]+}.phased.vcf'
	output: 'vcf/ashk.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.phased.vcf'
	log: 'vcf/ashk.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.phased.vcf.log'
	run:
		sample = role_to_sampleid[wildcards.individual]
		shell('vcf-subset -c {sample} {input} > {output} 2> {log}')


rule new_genetic_map:
	input: 'download/1000GP_Phase3/genetic_map_chr{chromosome,[0-9]+}_combined_b37.txt'
	output: 'genmap/scaled-{factor,[0-9]+}/genetic_map_chr{chromosome,[0-9]+}.txt'
	message: 'Multiplying genetic distances by {wildcards.factor}'
	shell: 'awk \'NR==1 {{print}} NR>1 {{print $1, $2*{wildcards.factor}, $3*{wildcards.factor} }}\' {input} > {output}'


## Create the VCF file for the "virtual child"

rule copy_sim_father_mother:
	input: 'vcf/ashk.{individual,(mother|father)}.chr{chromosome,[0-9]+}.phased.vcf'
	output: 'vcf/sim.{individual,(mother|father)}.chr{chromosome,[0-9]+}.phased.vcf'
	shell: 'ln -fsrv {input} {output}'


rule create_artificial_child:
	input:
		mothervcf='vcf/sim.mother.chr{chromosome,[0-9]+}.phased.vcf',
		fathervcf='vcf/sim.father.chr{chromosome,[0-9]+}.phased.vcf',
		genetic_map='genmap/scaled-10/genetic_map_chr{chromosome,[0-9]+}.txt'
	output:
		childvcf='vcf/sim.child.chr{chromosome,[0-9]+}.phased.vcf',
		motherrecomb='sim/mother.chr{chromosome,[0-9]+}.true.recomb',
		fatherrecomb='sim/father.chr{chromosome,[0-9]+}.true.recomb'
	log: 'vcf/sim.child.chr{chromosome,[0-9]+}.phased.vcf.log'
	message: 'Sampling artifical child'
	run: 
		sample = role_to_sampleid['child']
		shell('{artificial_child} {input.genetic_map} {input.mothervcf} {input.fathervcf} {sample} {output.motherrecomb} {output.fatherrecomb} > {output.childvcf} 2>{log}')


rule merge_artificial_trio:
	input:
		mothervcf='vcf/sim.mother.chr{chromosome,[0-9]+}.phased.vcf.gz',
		mothervcftbi='vcf/sim.mother.chr{chromosome,[0-9]+}.phased.vcf.gz.tbi',
		fathervcf='vcf/sim.father.chr{chromosome,[0-9]+}.phased.vcf.gz',
		fathervcftbi='vcf/sim.father.chr{chromosome,[0-9]+}.phased.vcf.gz.tbi',
		childvcf='vcf/sim.child.chr{chromosome,[0-9]+}.phased.vcf.gz',
		childvcftbi='vcf/sim.child.chr{chromosome,[0-9]+}.phased.vcf.gz.tbi',
	output:
		vcf='vcf/sim.trio.chr{chromosome,[0-9]+}.phased.vcf',
	log: 'vcf/sim.trio.chr{chromosome,[0-9]+}.phased.vcf.log'
	shell:
		"vcf-merge {input.mothervcf} {input.fathervcf} {input.childvcf} > {output.vcf} 2> {log}"


rule sim_fastas:
	input:
		ref=reference,
		vcf='vcf/sim.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.phased.vcf',
	output:
		hap1='sim/sim.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.haplotype1.true.fasta',
		hap2='sim/sim.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.haplotype2.true.fasta',
	log: 'sim/sim.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.true.log'
	message: 'Creating true haplotypes {output}'
	run:
		sample = role_to_sampleid[wildcards.individual]
		shell('mkdir -p sim/tmp')
		shell('{genomesimulator} -c {wildcards.chromosome} {input.vcf} {input.ref} sim/tmp > {log} 2>&1')
		shell('mv sim/tmp/{sample}.chr{wildcards.chromosome}.1.fasta {output.hap1}')
		shell('mv sim/tmp/{sample}.chr{wildcards.chromosome}.2.fasta {output.hap2}')


rule copy_gen_map_ashk:
	input: 'download/1000GP_Phase3/genetic_map_chr{chromosome,[0-9]+}_combined_b37.txt'
	output: 'genmap/ashk.chr{chromosome,[0-9]+}.map'
	shell: 'ln -fsrv {input} {output}'


rule copy_gen_map_sim:
	input: 'genmap/scaled-10/genetic_map_chr{chromosome,[0-9]+}.txt'
	output: 'genmap/sim.chr{chromosome,[0-9]+}.map'
	shell: 'ln -fsrv {input} {output}'


rule create_trio_ped:
	output: 'trio.ped'
	shell: 'echo family HG002 HG003 HG004 0 0 > {output}'


## Trio phasing rules are unused

rule trio_whatshap:
	input:
		ref=reference,
		vcf='vcf/{dataset,[a-z]+}.trio.chr{chromosome}.unphased.vcf',
		bamm='bam/{dataset,[a-z]+}.{platform,[a-z]+}.mother.chr{chromosome,[0-9]+}.cov{coverage,([0-9]+|all)}.bam',
		bamf='bam/{dataset,[a-z]+}.{platform,[a-z]+}.father.chr{chromosome,[0-9]+}.cov{coverage,([0-9]+|all)}.bam',
		bamc='bam/{dataset,[a-z]+}.{platform,[a-z]+}.child.chr{chromosome,[0-9]+}.cov{coverage,([0-9]+|all)}.bam',
		genmap='genmap/{dataset,[a-z]+}.chr{chromosome,[0-9]+}.map',
		ped='trio.ped'
	output:
		vcf='phased/whatshap/trio/{dataset,[a-z]+}.{platform,[a-z]+}.trio.chr{chromosome,[0-9]+}.cov{coverage,([0-9]+|all)}.vcf',
		#recomb='whatshap/trio/{dataset,[a-z]+}.{platform,[a-z]+}.chr{chromosome,[0-9]+}.cov{coverage,([0-9]+|all)}.recomb',
	log: 'phased/whatshap/trio/{dataset,[a-z]+}.{platform,[a-z]+}.chr{chromosome,[0-9]+}.cov{coverage,([0-9]+|all)}.log'
	# TODO
	# --recombination-list {output.recomb}
	shell:
		'{time} {whatshap} phase --reference {input.ref} --ped {input.ped} --genmap {input.genmap} --chromosome {wildcards.chromosome} -o {output.vcf} {input.vcf} {input.bamm} {input.bamf} {input.bamc} >& {log}'


rule split_trio_vcf:
	input:
		vcf='phased/whatshap/trio/{dataset,[a-z]+}.{platform,[a-z]+}.trio.chr{chromosome,[0-9]+}.cov{coverage,([0-9]+|all)}.vcf'
	output:
		vcf='phased/whatshap/trio/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.vcf'
	log: 'phased/whatshap/trio/{dataset}.{platform}.{individual}.chr{chromosome}.cov{coverage}.vcf.log'
	run:
		sample = role_to_sampleid[wildcards.individual]
		shell('vcf-subset -c {sample} {input} > {output} 2> {log}')


rule simulate_pacbio_reads:
	output:
		fastq='sim/sim.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.haplotype{hap,[12]}.fastq.gz',
		maf='sim/sim.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.haplotype{hap,[12]}.maf.gz'
	input:
		sample='fastq/ashk.pacbio.{individual}.chr{chromosome}.cov2.fastq',
		haplotype='sim/sim.{individual}.chr{chromosome}.haplotype{hap}.true.fasta'
	log: 'sim/sim.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.haplotype{hap,[12]}.fastq.log'
	run:
		coverage = 30
		halfcoverage = coverage / 2
		seed = abs(hash(output.fastq))
		shell('mkdir -p sim/tmp')
		shell('time (pbsim --seed {seed} --prefix sim/tmp/sim.{wildcards.individual}.chr{wildcards.chromosome}.haplotype{wildcards.hap} --depth {halfcoverage} --sample-fastq {input.sample} {input.haplotype}) > {log} 2>&1')
		shell('awk \'NR%4==1 {{printf("%s_HAP{wildcards.hap}\\n",$0)}} NR%4!=1 {{print}}\' sim/tmp/sim.{wildcards.individual}.chr{wildcards.chromosome}.haplotype{wildcards.hap}_0001.fastq | {gzip} > {output.fastq}')
		shell('cat sim/tmp/sim.{wildcards.individual}.chr{wildcards.chromosome}.haplotype{wildcards.hap}_0001.maf | {gzip} > {output.maf}')
		shell('rm -f sim/tmp/sim.{wildcards.individual}.chr{wildcards.chromosome}.haplotype{wildcards.hap}_*')


rule bwa_mem_single_end_pacbio:
	input: 
		fastq='sim/sim.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.haplotype{hap,[12]}.fastq.gz',
		ref=reference,
		indexed=reference + '.bwt'
	output: 'sim/sim.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.haplotype{hap,[12]}.bam'
	log: 'sim/sim.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.haplotype{hap,[12]}.bam.log'
	threads: 16
	run:
		sample = role_to_sampleid[wildcards.individual]
		shell('{time} bwa mem -x pacbio -t {threads} {input.ref} {input.fastq} | samtools view -u - | picard AddOrReplaceReadGroups {picard_tmpdir_switch} VALIDATION_STRINGENCY=LENIENT I=/dev/stdin O={output} ID={wildcards.individual}{wildcards.hap} LB=UNKNOWN PL=PACBIO PU=UNKNOWN SM={sample} 2> {log}')


rule merge_hap_bams:
	input:
		bam1='sim/sim.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.haplotype1.bam',
		bam2='sim/sim.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.haplotype2.bam'
	output: 
		bam='bam/sim.pacbio.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.covall.bam',
		bai='bam/sim.pacbio.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.covall.bai'
	log: 'bam/sim.pacbio.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.covall.bam.log'
	priority: 5
	message: 'Merging haplotype-specific BAM files to create {output.bam}'
	shell: '{time} picard -Xmx8g MergeSamFiles {picard_tmpdir_switch} VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=50000 SORT_ORDER=coordinate CREATE_INDEX=true CREATE_MD5_FILE=true I={input.bam1} I={input.bam2} O={output.bam} >& {log}'

################## we don't need until here.
## Rules for running the phasing tools.
## Some have more than one rule if pre- and/or postprocessing is needed

rule gatk_read_backed_phasing:
	input:
		bam='bam/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.bam',
		bai='bam/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.bai',
		vcf='vcf/{dataset,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.unphased.vcf',
		#ref=reference,
		#dictfile=reference.replace('.fasta', '.dict'),
		gatk_jar=gatk_jar
	output:
		vcf='phased/read-backed-phasing/noindels/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.vcf',
	log: 'phased/read-backed-phasing/noindels/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.log'
	message: 'Running GATK\'s ReadBackedPhasing on {input.bam}'
	threads: 1  # ReadBackedPhasing cannot run with more than one thread
	run:
		#input_line = ' '.join('-I ' + bam for bam in input.bams)
		ref = platform_to_ref[wildcards.platform]
		dictfile=platform_to_ref[wildcards.platform].replace('.fasta', '.dict')
		shell(r"""
		{time} java -jar {gatk_jar} -T ReadBackedPhasing \
			-R {ref} \
			-I {input.bam} \
			-L {input.vcf} \
			-V {input.vcf} \
			--phaseQualityThresh 1 --min_base_quality_score 1 \
			-o {output.vcf} >& {log}""")


# The hapcompass rule is unused because the tool runs out of memory (tried 200GB)
hapcompass_out = 'hapcompass/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}'

rule hapcompass:
	output:
		solution=hapcompass_out + '_MWER_solution.txt',
		frags=temp(hapcompass_out + '_frags.txt'),
		phsolution=temp(hapcompass_out + '_phasedSolution.txt'),
		reads_sam=temp(hapcompass_out + '_reads.sam'),
		rrsam=temp(hapcompass_out + '_reduced_representation.sam'),
		rrvcf=temp(hapcompass_out + '_reduced_representation.vcf'),
	params:
		prefix=hapcompass_out
	log:
		hapcompass_out + '.log'
	input:
		bam='bam/{dataset}.{platform}.{individual}.chr{chromosome}.cov{coverage}.bam',
		bai='bam/{dataset}.{platform}.{individual}.chr{chromosome}.cov{coverage}.bai',
		vcf='vcf/{dataset}.{individual}.chr{chromosome}.unphased.vcf',
	shell:
		"""{time} {hapcompass} --bam {input.bam} --vcf {input.vcf} -o {params.prefix} >& {log}"""


hapcut_out = '' + dataset_pattern + '.cov{coverage,([0-9]+|all)}'

rule extract_hairs:
	"""hapCUT-specific pre-processing"""
	output:
		txt='hairs/{indelsornot,(indels|noindels)}/' + hapcut_out + '.txt'
	input:
		#sref=reference,
		bam='bam/{dataset}.{platform}.{individual}.chr{chromosome}.cov{coverage}.bam',
		bai='bam/{dataset}.{platform}.{individual}.chr{chromosome}.cov{coverage}.bai',
		vcf='vcf/{dataset}.{individual}.chr{chromosome}.unphased.vcf',
	log:
		'hairs/{indelsornot,(indels|noindels)}/' + hapcut_out + '.log'
	run:
		extra = ' --indels 1' if wildcards.indelsornot == 'indels' else ''
		extra_quality = ' --noquality 10 ' if wildcards.platform == 'pacbio' else ''
		ref = platform_to_ref[wildcards.platform]
		shell("{time} {extract_hairs} --ref {ref} {extra} {extra_quality} --VCF {input.vcf} --bam {input.bam} --longreads 1 --maxIS 600 > {output.txt} 2> {log}")


rule hapcut:
	output:
		txt='phased/hapcut/{indelsornot,(indels|noindels)}/' + hapcut_out + '.txt'
	input:
		txt='hairs/{indelsornot}/' + hapcut_out + '.txt',
		vcf='vcf/{dataset}.{individual}.chr{chromosome}.unphased.vcf',
	log: 'phased/hapcut/{indelsornot}/' + hapcut_out + '.phase.log'
	shell:
		"{time} {hapcut} --fragments {input.txt} --VCF {input.vcf} --longreads 1 --output {output.txt} >& {log}"


# We use the extractHAIRS distributed with hapCUT as the one distributed with
# hapCUT 2 gives unusable output (perhaps the same problem as the one that
# we fixed in the hapCUT 1 version with SAM = and X operators)
rule hapcut2:
	output:
		txt='phased/hapcut2/{indelsornot,(indels|noindels)}/' + hapcut_out + '.txt'
	input:
		txt='hairs/{indelsornot}/' + hapcut_out + '.txt',
		vcf='vcf/{dataset}.{individual}.chr{chromosome}.unphased.vcf',
	log: 'phased/hapcut2/{indelsornot}/' + hapcut_out + '.phase.log'
	shell:
		"{time} {hapcut2} --fragments {input.txt} --VCF {input.vcf} --longreads 1 --output {output.txt} >& {log}"


rule hapcut_to_vcf:
	output:
		vcf='phased/hapcut/{indelsornot,(indels|noindels)}/' + hapcut_out + '.vcf'
	input:
		vcf='vcf/{dataset}.{individual}.chr{chromosome}.unphased.vcf',
		hapcut='phased/hapcut/{indelsornot}/' + hapcut_out + '.txt'
	log: 'phased/hapcut/{indelsornot}/' + hapcut_out + '.vcf.log'
	shell:
		"{time} {whatshap} hapcut2vcf {input.vcf} {input.hapcut} > {output.vcf} 2> {log}"


rule hapcut2_to_vcf:
	output:
		vcf='phased/hapcut2/{indelsornot,(indels|noindels)}/' + hapcut_out + '.vcf'
	input:
		vcf='vcf/{dataset}.{individual}.chr{chromosome}.unphased.vcf',
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
		bam='bam/{dataset}.{platform}.{individual}.chr{chromosome}.cov{coverage}.bam',
		bai='bam/{dataset}.{platform}.{individual}.chr{chromosome}.cov{coverage}.bai',
		vcf='vcf/{dataset}.{individual}.chr{chromosome}.unphased.vcf.gz',
		tbi='vcf/{dataset}.{individual}.chr{chromosome}.unphased.vcf.gz.tbi',
	log:
		'phased/phaser/{indelsornot}/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.log'
	run:
		extra = ' --include_indels 1' if wildcards.indelsornot == 'indels' else ' --include_indels 0'
		sample = role_to_sampleid[wildcards.individual]
		shell("{time} {phaser}{extra} --bam {input.bam} --write_vcf 1 --gw_phase_vcf 2 --pass_only 0 --vcf {input.vcf} --sample {sample} --mapq 1 --baseq 1 --paired_end 0 --o {params.base} >& {log}")
		shell("gunzip -f {params.base}.vcf.gz")
		#shell(r"sed -i.orig '/^#/!s|\bGT:PG:|PG:GT:|;/^#/!s|:PI:|:PS:|' {params.base}.vcf")


rule whatshap_norealign:
	input:
		bam='bam/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.bam',
		vcf='vcf/{dataset,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome,[0-9]+}.unphased.vcf'
	output: 'phased/whatshap-norealign/noindels/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.vcf'
	log: 'phased/whatshap-norealign/noindels/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.log'
	shell: '{time} {whatshap} phase {input.vcf} {input.bam} > {output} 2> {log}'


rule whatshap_noindels:  # with re-alignment
	input:
		#ref=platform_to_ref[platform],
		bam='bam/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.bam',
		vcf='vcf/{dataset,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome, (all|[0-9]+)}.unphased.vcf'
	output:
		vcf='phased/whatshap/noindels/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.vcf',
		log='phased/whatshap/noindels/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.log'
	run: 
		ref = platform_to_ref[wildcards.platform]
		shell('{time} {whatshap} phase --reference {ref} {input.vcf} {input.bam} > {output.vcf} 2> {output.log}')


rule whatshap_indels:  # with re-alignment
	input:
		#ref=platform_to_ref[platform],
		bam='bam/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.bam',
		vcf='vcf/{dataset,[a-z]+}.{individual,(mother|father|child)}.chr{chromosome, (all|[0-9]+)}.unphased.vcf'
	output:
		vcf='phased/whatshap/indels/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.vcf',
		log='phased/whatshap/indels/' + dataset_pattern + '.cov{coverage,([0-9]+|all)}.log'
	run:
		ref = platform_to_ref[wildcards.platform]
		shell('{time} {whatshap} phase --indels --reference {ref} {input.vcf} {input.bam} > {output.vcf} 2> {output.log}')


## Evaluation: compare phasing results to ground truth, get VCF statistics

#TODO: did not update evaluation and further rules.
rule evaluate_phasing_tool:
	input:
		truth='vcf/{dataset}.{individual}.chr{chromosome}.phased.vcf',
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
			print('dataset', 'individual', 'indels', 'coverage', 'components', sep='\t', file=f)
			indels = int(wildcards.indelsornot == 'indels')
			print(wildcards.dataset, wildcards.individual, indels, wildcards.coverage, components, sep='\t', file=f)



## Create two tables with evaluation results

def all_eval_paths(extension):
	l1 = expand('eval/{program}/{dataset}.pacbio.{individual}.chr1.cov{coverage}' + extension ,
		program=['hapcut/indels', 'hapcut/noindels', 'hapcut2/indels', 'hapcut2/noindels',
			'read-backed-phasing/noindels',
			'whatshap-norealign/noindels', 'whatshap/noindels', 'whatshap/indels', 'whatshap/trio'],
		dataset=datasets, individual=individuals, coverage=coverage),
	l2 = expand('eval/phaser/{indels}/{dataset}.pacbio.{individual}.chr1.cov{coverage}' + extension,
		indels=['indels', 'noindels'], dataset=datasets, individual=individuals, coverage=coverage)
	# Exclude because they did not finish within 24h
	l2 = tuple(s for s in l2 if not ('/ashk.' in s and '.covall.' in s))
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
		components=expand('eval/components/{indels}/{dataset}.pacbio.{individual}.chr1.cov{coverage}.components',
			indels=['indels', 'noindels'], dataset=datasets, individual=individuals, coverage=coverage)
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


rule index_reference:
	output:
		'{path}.fasta.bwt'
	input:
		'{path}.fasta'
	shell:
		"bwa index {input}"


rule CreateSequenceDict:
	output: '{base}.dict'
	input: '{base}.fasta'
	shell:
		"picard CreateSequenceDictionary R={input} O={output}"
