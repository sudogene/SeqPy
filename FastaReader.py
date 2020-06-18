from Seq import *

class FastaReader:
    _classes = {'Seq': Seq, 'DNA': DNA, 'RNA': RNA, 'Protein': Protein}

    def __init__(self, filename, seq_type='Seq'):
        if not filename.endswith('.fasta') or seq_type not in FastaReader._classes.keys():
            raise ValueError("Invalid input")

        self._sequences = []
        with open(filename, 'r') as f:
            line = f.readline()
            seq = ''
            header = ''
            while line:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if seq:
                        self._sequences.append((header, FastaReader._classes[seq_type](seq)))
                        seq = ''
                        header = line
                    else:
                        header = line
                else:
                    seq += line
                line = f.readline()
            self._sequences.append((header, FastaReader._classes[seq_type](seq)))

    def __repr__(self):
        return repr(self._sequences)

    def __len__(self):
        return len(self._sequences)
    
    def __getitem__(self, other):
        return self._sequences[other]


class FastaObject:
    def __init__(self, string):
        