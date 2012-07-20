import sys, os
from argparse import ArgumentParser

def merge(left_file, right_file):
    try:
        left = open(left_file, 'r')
        right = open(right_file, 'r')
        filebase, ext = os.path.splitext(left_file)
        combined_file = filebase.split('_')[0]
        combined = open(combined_file, 'w')
        counter = 0
        seqs = []
        for line in left:
            if '@' in line:
                combined.write('@' + str(counter) + '\n')
                r = right.readline()
            else:
                combined.write(line)
                r = right.readline()
                counter += 1
                combined.write('@' + str(counter) + '\n')
                combined.write(r)
                counter += 1
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
    merge(args.read_file_1, args.read_file_2)
    
if __name__ == '__main__':
    sys.exit(main())
