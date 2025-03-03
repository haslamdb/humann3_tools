

# TODO: Issues and Fixes for humann3_tools

## 🛠️ Bugs to Fix
- [ ] --paired flag not working because kneaddata fails with paried
- [ ] Move kneaddata output files and humann output files to new directories then get rid of intermediate files
- [ ] Change --paired flag to --un (unpaired)
- [ ] Add flag for --decontaminate-pairs strict

`

## 🔧 Features to Improve
- [ ] Check if intermediate files already exist before performing steps in the analysis
        e.g. if Kneaddata output files already exist, skip kneaddata


```bash

# no paired flag - as it interferes with kneaddata
	
humann3-tools --log-file test.log \
	--log-level DEBUG \
	--output-dir test.out \
	--output-prefix test- \
	--run-preprocessing \
	--kneaddata-dbs  /home/david/Databases/KneadDataDB /home/david/Databases/BT2ContaminantDB \
	--threads 24 \
	--seq-dir /home/david/Data/TrimmedMSSFiles/ \
	--input-fastq /home/david/Data/TrimmedMSSFiles/D19G_R1.fastq /home/david/Data/TrimmedMSSFiles/D19G_R2.fastq \
	--sample-key TestKey.csv \
	--samples-file ExampleFileList.txt \
	--use-metadata \
	--pathway-dir test.out \
	--gene-dir test.out \
	--annotations-dir annotations \
	--skip-downstream
