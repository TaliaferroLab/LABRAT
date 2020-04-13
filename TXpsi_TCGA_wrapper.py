#from subprocess import Popen, PIPE, STDOUT
import subprocess

cancers = ['BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'ESCA', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'PAAD', 'PCPG', 'PRAD', 'READ', 'STAD', 'THCA', 'THYM', 'UCEC']

gff = '/beevol/home/taliaferro/Annotations/hg38/Gencode28/gencodecomprehensive.v28.gff3'

for idx, cancer in enumerate(cancers):
	command = ['python', '/beevol/home/taliaferro/Scripts/TXpsi_TCGA.py', '--gff', gff, '--cancerdir', cancer]
	print 'Analyzing {0}, sample {1} of {2}...'.format(cancer, idx + 1, len(cancers))
	subprocess.call(command)
