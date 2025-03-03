humann3-tools --log-file test.log \
	--log-level DEBUG \
	--output-dir test.out \
	--output-prefix test- \
	--run-preprocessing \
	--paired \
	--kneaddata-dbs  /home/david/Databases/KneadDataDB /home/david/Databases/BT2ContaminantDB \
	--threads 12 \
  	--kneaddata-output-dir /path/to/kneaddata/results \
  	--humann3-output-dir /path/to/humann3/results \
	--seq-dir /home/david/Data/TrimmedMSSFiles/ \
	--input-fastq /home/david/Data/TrimmedMSSFiles/D19G_R1.fastq /home/david/Data/TrimmedMSSFiles/D19G_R2.fastq \
	--sample-key TestKey.csv \
	--samples-file ExampleFileList.txt \
	--use-metadata \
	--pathway-dir test.out \
	--gene-dir test.out \
	--annotations-dir annotations \
	--skip-downstream
	

# no paired flag - as it interferes with kneaddata
	
humann3-tools --log-file test.log \
	--log-level DEBUG \
	--output-dir test.out \
	--output-prefix test- \
	--run-preprocessing \
	--kneaddata-dbs  /home/david/Databases/KneadDataDB /home/david/Databases/BT2ContaminantDB \
	--threads 24 \
  	--kneaddata-output-dir ~/Data/TrimmedMSSFiles \
  	--humann3-output-dir ~/Documents/Alignments/Humann3Alignments \
	--sample-key TestKey.csv \
	--use-metadata \
	--seq-dir /home/david/Data/TrimmedMSSFiles/ \
	--pathway-dir test.out \
	--gene-dir test.out \
	--group-col "Group" \
	--annotations-dir annotations 
	
	--samples-file ExampleFileList.txt \
	--input-fastq /home/david/Data/TrimmedMSSFiles/D19G_R1.fastq /home/david/Data/TrimmedMSSFiles/D19G_R2.fastq \


