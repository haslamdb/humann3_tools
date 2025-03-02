

# TODO: Issues and Fixes for MyPackage

## 🛠️ Bugs to Fix
- [ ] --paired flag not working in from cli because kneaddata doesn't acccept paired flag
`

## 🔧 Features to Improve
- [ ] There are redundancies in the sample input. For example, below is a command
        that contains several ways to get the input files, but I think it's just using the --input-fastq flag.
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
