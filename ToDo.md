

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


	
humann3-tools --log-file test.log \
	--log-level DEBUG \
	--output-dir test.out2 \
	--output-prefix test- \
	--run-preprocessing \
	--kneaddata-dbs  /home/david/Databases/KneadDataDB /home/david/Databases/BT2ContaminantDB \
	--threads 24 \
	--sample-key TestKey.csv \
	--use-metadata \
	--seq-dir /home/david/Data/TrimmedMSSFiles/ \
  	--kneaddata-output-dir ~/Data/TrimmedMSSFiles \
  	--humann3-output-dir ~/Documents/Alignments/Humann3Alignments \
	--pathway-dir test.out2 \
	--gene-dir test.out2 \
	--group-col "Group" \
	--paired \
	--decontaminate-pairs strict \
	--annotations-dir annotations 
	
