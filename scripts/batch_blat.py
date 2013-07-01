from zoom import *
import glob

QUERY_DIR = '/home/ariyaratnep/shaojiang/peta/rnaseq/hg19/SRX011545/'

def run(chr_dir):
	for f in glob.glob(chr_dir + "/chr*"):
		cmd = '/home/ariyaratnep/shaojiang/blat/blat -noHead %s %s/both.fa -ooc=/home/ariyaratnep/shaojiang/blat/11.ooc %s.psl' % (f, QUERY_DIR, f)
		print cmd
		runInShell(cmd)
	files = ''
	for f in glob.glob(chr_dir + "/chr*.psl"):
		files += f + ' '
	cmd = 'cat %s > %s/both.genome.psl' % (files, QUERY_DIR)
	print cmd
	runInShell(cmd)
run('/home/ariyaratnep/shaojiang/peta/rnaseq/hg19/')
