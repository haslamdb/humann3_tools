

# TODO: Issues and Fixes for humann3_tools

## 🛠️ Bugs to Fix
- [ ] move kneaddata and humann3 output files to a different directory then delete the temp directory
- [ ] in cli.py , if args.use_metadata is True but no --seq-dir is provided, the current code does not execute metadata-based file collection.

`

## 🔧 Features to Improve
- [ ] Check if intermediate files already exist before performing steps in the analysis
        e.g. if Kneaddata output files already exist, skip kneaddata
- [ ] Add genedirectory, pathwaydirectory, etc flags to options section in README.md
- [ ] incorporate into the README.md new join_unstratify_humann_output command


```bash


	
# not parallel	
	humann3-tools --run-preprocessing \
	--log-file test.log \
	--log-level DEBUG \
	--output-dir test.out2 \
	--output-prefix test- \
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


# parallel
humann3-tools 	--run-preprocessing \
	--log-file test.log \
	--log-level DEBUG \
	--output-dir test.out2 \
	--output-prefix test- \
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
	--use-parallel \
	--threads-per-sample 4 \
	--decontaminate-pairs strict \
	--annotations-dir annotations 

# join, normalize, split stratify only
humann3-tools --join-only \
  --sample-key TestKey.csv \
  --pathway-dir /path/to/pathways \
  --gene-dir /path/to/genes \
  --output-dir joined_output \
  --units cpm \
  --no-interactive
