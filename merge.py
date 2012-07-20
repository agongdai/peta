import sys, os
from argparse import ArgumentParser

def merge_fq(left_file, right_file):
    try:
        left = open(left_file, 'r')
        right = open(right_file, 'r')
        filebase, ext = os.path.splitext(left_file)
        combined_file = filebase.split('_')[0] + '.fastq'
        combined = open(combined_file, 'w')
        counter = 0
        seqs = []
        line = left.readline()
        while line:
            if '@' in line:
#                print line
                combined.write('@' + str(counter) + '\n')
                l = left.readline()
                combined.write(l)
                l = left.readline()
                combined.write(l)
                l = left.readline()
                combined.write(l)
                counter += 1
                
                r = right.readline()
                combined.write('@' + str(counter) + '\n')
                r = right.readline()
                combined.write(r)
                r = right.readline()
                combined.write(r)
                r = right.readline()
                combined.write(r)
                counter += 1
            line = left.readline()
                
        left.close()
        right.close()
        combined.close()
    except:
        print("IOError")

def main():
    parser = ArgumentParser()
    parser.add_argument('read_file_1', help='Left mates file')
    parser.add_argument('read_file_2', help='Right mates file')
    args = parser.parse_args()
    merge_fq(args.read_file_1, args.read_file_2)
    
if __name__ == '__main__':
    sys.exit(main())
