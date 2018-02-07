#from subprocess import Popen, PIPE, STDOUT
import subprocess

cancers = ['BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'ESCA', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'PAAD', 'PCPG', 'PRAD', 'READ', 'STAD', 'THCA', 'THYM', 'UCEC']

gff = '/vol3/home/taliaferro/Annotations/hg38/gencode.hg38comprehensive.gff3'

for idx, cancer in enumerate(cancers):
	command = ['python', '/vol3/home/taliaferro/Scripts/TXpsi_TCGA.py', '--gff', gff, '--cancerdir', cancer]
	print 'Analyzing {0}, sample {1} of {2}...'.format(cancer, idx + 1, len(cancers))
	subprocess.call(command)
