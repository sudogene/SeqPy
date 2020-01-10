class Seq:
    def __init__(self, seq):
        self._seq = seq


    def __repr__(self):
        if len(self._seq) > 21:
            return f'{self.__class__.__name__}({self._seq[:21]}...{self._seq[-3:]})'
        return f'{self.__class__.__name__}({self._seq})'


    def __str__(self):
        return self._seq


    def __add__(self, other):
        return self.__class__(self._seq + other._seq)


    def __mul__(self, other):
        return self.__class__(self._seq * other)


    def __contains__(self, item):
        return item in self._seq


    def __len__(self):
        return len(self._seq)


    def __getitem__(self, key):
        return self.__class__(self._seq[key])


    def __eq__(self, other):
        return self._seq == other._seq


    def _class_check(self, other):
        '''
        Ensures other class is a valid class for subsequent computation
        '''
        return isinstance(other, self.__class__) or \
                issubclass(other.__class__, self.__class__) or \
                issubclass(self.__class__, other.__class__) or \
                type(other) == str


    def count(self, pattern):
        if not self._class_check(pattern):
            raise TypeError("Invalid input type")
        return self._seq.count(str(pattern))


    def count_overlap(self, pattern):
        if not self._class_check(pattern):
            raise TypeError("Invalid input type")
        
        pattern = str(pattern)
        count = start = 0
        while True:
            start = self._seq.find(pattern, start) + 1
            if start > 0:
                count += 1
            else:
                return count


class DNA(Seq):
    _transcribe = {'A':'U', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
    _complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
    _valid = {'A', 'T', 'G', 'C', 'N'}

    def __init__(self, seq): 
        for base in seq:
            if base not in self._valid:
                raise ValueError("Invalid sequence")       
        super().__init__(seq)


    def transcribe(self):
        seq = ''.join([self._transcribe[b] for b in self._seq])
        return RNA(seq)


    def t(self):
        return self.transcribe()


    def complement(self):
        seq = ''.join([self._complement[b] for b in self._seq])
        return DNA(seq)


    def c(self):
        return self.complement()


    def reverse_complement(self):
        return self.complement()[::-1]


    def rc(self):
        return self.reverse_complement()


    def gc(self):
        val = (self._seq.count('G') + self._seq.count('C')) / len(self._seq)
        return round(val, 2)
    

    def codons(self):
        return tuple(map(''.join, zip(*[iter(self._seq)]*3)))
    

    def frames(self, complete=True):
        result = {}
        indices = (1, 2, 3, -1, -2, -3) if complete else (1, 2, 3)
        for frame in indices:
                if frame < 0:
                    result[frame] = self.c()[((frame * -1) - 1):]
                else:
                    result[frame] = self[(frame - 1):]
        return result


class RNA(Seq):
    _valid = {'A', 'U', 'G', 'C', 'N'}

    def __init__(self, seq):
        for base in seq:
            if base not in self._valid:
                raise ValueError("Invalid sequence")
        super().__init__(seq)
