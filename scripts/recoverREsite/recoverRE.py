import sys
import argparse

def meanPhred(scores):
    phredScores = {'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, "'": 6, '(': 7, ')': 8,
                '*': 9, '+': 10, ',': 11, '-': 12, '.': 13, '/': 14, '0': 15, '1': 16,
                '2': 17, '3': 18, '4': 19, '5': 20, '6': 21, '7': 22, '8': 23, '9': 24,
                ':': 25, ';': 26, '<': 27, '=': 28, '>': 29, '?': 30, '@': 31, 'A': 32,
                'B': 33, 'C': 34, 'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39, 'I': 40,
                'J': 41, 'K': 42}
    inv_dict = {v: k for k, v in phredScores.items()}

    summ = 0.0
    for char in scores:
        summ += phredScores[char]
    mean = int(summ/len(scores))
    mean_enc = inv_dict[mean]
    
    return(mean_enc * 5)

def parseFile(inputname, start, end, truesite, outname, fixscore):
    resites = []
    with open(inputname) as fastafile, open(outname, 'w') as outfile:
        onoff = False
        onoff_score = False
        for i, line in enumerate(fastafile):
            l = line.replace('\n', '')
            newseq = ''
            if onoff == True:
                if len(l) > 10:
                    if start == 1: #Read 2
                        resite = l[:end]
                        therest = l[end:]
                        newseq = truesite+therest
                    elif start != 1:
                        resite = l[start-1:end]
                        barcode = l[0:start-1]
                        therest = l[end:]
                        newseq = barcode+truesite+therest
                else:
                    newseq = l
                assert len(newseq) == len(l) 
                outfile.write(newseq+'\n')
                onoff = False

            elif onoff_score == True:
                if len(l) > 10:
                    if start == 1:
                        resite_phred = l[:end]
                        assert len(resite_phred) == 5
                        therest_phred = l[end:]
                        newScore = meanPhred(l)
                        newseq = newScore+therest_phred
                
                    elif start != 1:
                        resite_phred = l[start-1:end]
                        barcode_phred = l[0:start-1]
                        therest_phred = l[end:]
                        newScore = meanPhred(l)
                        newseq = barcode_phred+newScore+therest_phred
                else:
                    newseq = l
                assert len(newseq) == len(l)
                outfile.write(newseq + '\n')
                onoff_score = False
            
            elif "@NB" in l:
                outfile.write(l + '\n')
                onoff = True
            elif l == "+":
                outfile.write(l + '\n')
                onoff_score = True
            else:
                outfile.write(l + '\n')

def main():
    parser = argparse.ArgumentParser(description = "Recover(fix) RE sites in ddRAD data to facilitate mapping and \
            SNP calling with ipyrad and Stacks")
    parser.add_argument("input", help="Location of read file")
    parser.add_argument("--truesite", help="True sequence of RE site")
    parser.add_argument("--startpos", help="Start position of RE site in read", default=6, type=int)
    parser.add_argument("--endpos", help="End position of RE site in read", default=10, type=int)
    parser.add_argument("--out", help="Location+name to output fixed file", default='./')
    parser.add_argument("--fixscore", help="Boolean. If '1' is supplied, the phred quality scores at the RE site will be adjusted to the mean quality of the read", default=0)

    args = parser.parse_args()
    
    if args.fixscore == '1':
        fixscore = True
    else:
        fixscore = False

    parseFile(args.input, int(args.startpos), int(args.endpos), args.truesite, args.out, fixscore)

if __name__ == '__main__':
    main()
