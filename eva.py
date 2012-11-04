import sys, os
from argparse import ArgumentParser

def differ(args):
	files = args.inputs.split(',')
	out_file = open(args.output, 'w')
	genes = {}
	for f in files:
		fp = open(f, 'r')
		for line in fp:
			line = line.strip()
			if line in genes:
				genes[line] += ',' + f
			else:
				genes[line] = f
		fp.close()
	categories = {}
	for key, value in genes.iteritems():
		out_file.write(key + ': ' + value  + '\n')
		if value in categories:
			categories[value] += 1
		else:
			categories[value] = 0
	for key, value in categories.iteritems():
		print key, value
	out_file.close()

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub command help')
    parser_differ = subparsers.add_parser('differ', help='differ multiple files')
    parser_differ.set_defaults(func=differ)
    parser_differ.add_argument('-i', required=True, help='the files to compare, seperated by a ","', dest='inputs')
    parser_differ.add_argument('-o', required=True, help='result file', dest='output')

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    sys.exit(main())
