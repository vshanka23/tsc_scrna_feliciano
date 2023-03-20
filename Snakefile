configfile: "cellranger_counts.yaml"

SAMPLES=sorted(glob_wildcards(config["SOURCE"]+"/{samples}"+config["EXT"]).samples)
SAMPLES_NF=[]
import re
for item in SAMPLES:
	var1=re.sub('_S[0-9]_L00[1-2]_[RI][1-2]','',str(item))
	SAMPLES_NF.append(var1)

SAMPLES_SET = set(SAMPLES_NF)
SAMPLES_ID = list(SAMPLES_SET)

rule all:
	input:
		expand(config["DEST"]+"/{SAMPLES_ID}.txt",SAMPLES_ID=SAMPLES_ID)

rule counting:
	output:
		config["DEST"]+"/{SAMPLES_ID}.txt"
	params:
		TRANSCRIPTOME=config["INDEX_REF"],
		SAMPLE="{SAMPLES_ID}",
		FASTQS=config["SOURCE"],
		OUTPUT_DIR=config["DEST"]
	resources: cpus=5, mem_mb=20000, time_min=1440
	shell:
		"""
		ml cellranger/7.0.1
		cd {params.OUTPUT_DIR}
		touch {params.SAMPLE}.txt
		cellranger count --id={params.SAMPLE} \
		--fastqs={params.FASTQS} \
		--sample={params.SAMPLE} \
		--transcriptome={params.TRANSCRIPTOME} \
		--force-cells 60000 \
		--localcores=5 \
		--localmem=20 \
		"""
