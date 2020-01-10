from Seq import *

class FastaReader:
    def __init__(self, filename, seq_type='Seq'):
        if not filename.endswith('.fasta') or seq_type not in ['Seq', 'DNA', 'RNA']:
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
                        if seq_type == 'Seq':
                            self._sequences.append((header, Seq(seq)))
                        elif seq_type == 'DNA':
                            self._sequences.append((header, DNA(seq)))
                        elif seq_type == 'RNA':
                            self._sequences.append((header, RNA(seq)))
                        seq = ''
                        header = line
                    else:
                        header = line
                else:
                    seq += line
                line = f.readline()
            if seq_type == 'Seq':
                self._sequences.append((header, Seq(seq)))
            elif seq_type == 'DNA':
                self._sequences.append((header, DNA(seq)))
            elif seq_type == 'RNA':
                self._sequences.append((header, RNA(seq)))

    def __repr__(self):
        return repr(self._sequences)

    def __len__(self):
        return len(self._sequences)
    
    def __getitem__(self, other):
        return self._sequences[other]
