from Seq import *
import re


class FastaReader:
    """
    FastaReader parses fasta-formatted files to store all 
    fasta records as FastaRecords. Can be iterated over in a loop
    to retrieve all records.
    """
    _extensions = {'txt', 'fasta', 'fa'}
    _classes = {'Seq': Seq, 'DNA': DNA, 'RNA': RNA, 'Protein': Protein}

    def __init__(self, filename, seq_type='Seq'):
        if filename.split('.')[-1] not in self._extensions or \
                seq_type not in self._classes.keys():
            raise ValueError("invalid input")

        self._sequences = []
        with open(filename, 'r') as f:
            line = f.readline()
            seq = ''
            header = ''
            while line:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if seq:
                        fo = FastaRecord(header, self._classes[seq_type](seq.upper()))

                        self._sequences.append(fo)
                        seq = ''
                        header = line
                    else:
                        header = line
                else:
                    seq += line
                line = f.readline()
            fo = FastaRecord(header, self._classes[seq_type](seq.upper()))
            self._sequences.append(fo)

    def __repr__(self):
        return f"FastaReader[{len(self._sequences)}]"

    def __len__(self):
        return len(self._sequences)

    def __getitem__(self, other):
        return self._sequences[other]


class FastaRecord:
    """
    FastaRecord contains a dictionary of the description of the 
    fasta entity, usually found within the header, and the sequence.
    
    The record can be queried like a dictionary, if needed, by following
    the standard dictionary functions such as keys(), values(), items(),
    and the dunder __getitem__ method. Sequence can be retrieved via seq.
    """

    def __init__(self, header, seq):
        self.desc = {'id': header.split()[0]}
        descs = list(re.findall(r"[^[]*\[([^]]*)]", header))
        for d in descs:
            k, v = d.split('=')
            self.desc[k] = v
        self.seq = seq

    def __repr__(self):
        return f"Fasta({self.desc['id']})"

    def __len__(self):
        return len(self.seq)

    def describe(self):
        d = []
        for k, v in self.desc.items():
            d.append(f"{k}={v}")
        print("\n".join(d))

    def __getitem__(self, k):
        try:
            return self.desc[k]
        except KeyError:
            raise KeyError("Call keys() for valid keys")

    def keys(self):
        return self.desc.keys()

    def values(self):
        return self.desc.values()

    def items(self):
        return self.desc.items()


if __name__ == '__main__':
    fr = FastaReader('sequence.txt', 'DNA')
