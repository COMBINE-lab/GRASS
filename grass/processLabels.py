import os
from subprocess import call
import pandas as pd

#run BLAST and return list of label files in order [TS->ASdb, AS->TSdb] 
def runBLAST(faFile1, faFile2, outdir):
	TSdb = os.path.sep.join([outdir, "TS", "db"])
	ASdb = os.path.sep.join([outdir, "AS", "db"])
	call(["makeblastdb", "-in", faFile1, "-dbtype", "nucl", "-out", TSdb])
	call(["makeblastdb", "-in", faFile2, "-dbtype", "nucl", "-out", ASdb])

	labelFile1 = os.path.sep.join([outdir, "TS.ASdb" ]) 
	labelFile2 = os.path.sep.join([outdir, "AS.TSdb" ]) 
	#run blastn for TS against AS
	call(["blastn", "-db", ASdb, "-query", faFile1, "-outfmt", "6", "-out", labelFile1, "-num_threads", "8"])
	call(["blastn", "-db", TSdb, "-query", faFile1, "-outfmt", "6", "-out", labelFile2, "-num_threads", "8"])

	return [labelFile1, labelFile2]

#read in the BLAST results and process them to get two way best match. This creates the seed file for Junto.
#Label files should be in order [TS->ASdb, AS->TSdb] 
def genFinalLabels(keys, labelFiles, finalLabelFile, outdir):
	blastOut = open(keys["seed_file"], 'w')
	if (len(labelFiles)==1):
		with open(labelFiles[0].strip(), 'r') as f:
			data = pd.read_table(f, header=None, usecols=[0,1], names=['query', 'subject'])
			table1 = data.set_index("query").to_dict()['subject']
		try:
			for key,value in table1.items():
				probability=1.0
				blastOut.write(key + "\t" + value + "\t" +str(probability)+ "\n")
		except KeyError:
			pass
	else:
		labelFiles[0] = labelFiles[0].strip()
		labelFiles[1] = labelFiles[1].strip()
		outfile1 = outdir + "best" + os.path.basename(labelFiles[0])
		outfile2 = outdir + "best" + os.path.basename(labelFiles[1])
		command = "cat " + labelFiles[0]  + " |  sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > " + outfile1
		call(command, shell=True)
		command = "cat " + labelFiles[1] + " |  sort -k2,2 -k12,12nr -k11,11n | sort -u -k2,2 --merge > " + outfile2
		call(command, shell=True)
		with open(outfile1, 'r') as f:
				data = pd.read_table(f, header=None, usecols=[0,1], names=['query', 'subject'])
				table1 = data.set_index("query").to_dict()['subject']
	
		with open(outfile2, 'r') as f:
				data = pd.read_table(f, header=None, usecols=[0,1], names=['query', 'subject'])
				table2 = data.set_index("subject").to_dict()['query']
		call(["rm", outfile1])
		call(["rm", outfile2])

		try:
			for key,value in table1.items():
				if key in table2:
					probability=1.0
					blastOut.write(key + "\t" + value + "\t" +str(probability)+ "\n")
		except KeyError:
				pass
	blastOut.close()
	call(["cp", keys["seed_file"], finalLabelFile])