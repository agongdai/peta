import sys

def merge(left_file, right_file):
    try:
        left = open(left_file, 'r')
        right = open(right_file, 'r')
        combined = open('reads.fa', 'w')
        counter = 0
        seqs = []
        for line in left:
            if '>' in line:
                combined.write('>' + str(counter) + '\n')
                r = right.readline()
            else:
                combined.write(line)
                r = right.readline()
                counter += 1
                combined.write('>' + str(counter) + '\n')
                combined.write(r)
                counter += 1
        left.close()
        right.close()
        combined.close()
    except:
        print("IOError")

def main():
    merge('/home/carl/Projects/trinityrnaseq_r2012-01-25/SRR097897_1.fasta','/home/carl/Projects/trinityrnaseq_r2012-01-25/SRR097897_2.fasta')
    
if __name__ == '__main__':
    sys.exit(main())
