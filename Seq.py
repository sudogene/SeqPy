import Seq_data

class Seq:
    '''
    Seq (Sequence) is an abstract Python string wrapper class with added
    functionality for biological sequences such as nucleotides and peptides.
    '''

    def __init__(self, seq):
        self._seq = seq.upper()


    def __repr__(self):
        LIMIT = 21
        if len(self._seq) > LIMIT:
            return f'{self.__class__.__name__}({self._seq[:LIMIT]}...{self._seq[-3:]})'
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
        if self.__class__ == other.__class__:
            return self._seq == other._seq
        else:
            return False


    def _class_check(self, other):
        '''
        Ensures other class is a valid class for subsequent computation
        '''
        return isinstance(other, self.__class__) or \
                issubclass(other.__class__, self.__class__) or \
                issubclass(self.__class__, other.__class__) or \
                type(other) == str


    def _map(self, mapper, seqclass=None):
        '''
        The generic map function for all sequence mapping, allows for a robust abstract method
        for all variations of mapping (mapper), and the return sequence subclass (seqclass).

        Examples of mapping includes complement, transcription, and translation.
        '''
        if not seqclass:
            seqclass = self.__class__
        try:
            new_seq = ''.join([mapper[base] for base in self._seq])
            return seqclass(new_seq)
        except KeyError:
            raise KeyError('The provided mapper does not contain some bases in the sequence.')
    

    def reverse(self):
        return self.__class__(self._seq[::-1])
    

    def count(self, pattern):
        if not self._class_check(pattern):
            raise TypeError('Invalid input type')
        return self._seq.count(str(pattern))


    def count_overlap(self, pattern):
        if not self._class_check(pattern):
            raise TypeError('Invalid input type')
        
        pattern = str(pattern)
        count = start = 0
        while True:
            start = self._seq.find(pattern, start) + 1
            if start > 0:
                count += 1
            else:
                return count


class AbstractNucleotide(Seq):
    '''
    AbstractNucleotide is the abstract class for nucleotide sequences including DNA and RNA.
    It has parent functions for finding complements, generating reading frames and codons.
    '''

    def __init__(self, seq, valid_checker):
        '''
        Each nucleotide type has a set of valid bases in the valid_checker argument.
        Most subclasses of AbstractNucleotide have their own default sets;
        e.g. {A, T, G, C, N} for DNA and {A, U, G, C, N} for RNA.

        Users can also insert their own set if desired.
  
        '''
        for base in seq:
            if base not in valid_checker:
                raise ValueError('Invalid sequence')       
        super().__init__(seq)


    def complement(self, mapper):
        seq = ''.join([mapper[b] for b in self._seq])
        return self.__class__(seq)


    def c(self):
        return self.complement()


    def gc(self):
        '''
        Returns GC content of the nucleotide sequence, rounded to 2 decimal places.
        
        >>> DNA("AATTGGCC").gc()
        0.67
        
        '''

        val = (super().count('G') + super().count('C')) / len(self._seq)
        return round(val, 2)
    

    def codons(self, keep_class=False):
        '''
        A codon is a nucleotide with 3 bases, read during translation event.
        Returns a tuple of codon strings by default, else a tuple of codon nucleotides.

        >>> DNA("AGGCCTGGGGTTTAATTAGCGTAGCGTTTACGATAT").codons()
        ('AGG', 'CCT', 'GGG', 'GTT', 'TAA', 'TTA', 'GCG', 'TAG', 'CGT', 'TTA', 'CGA', 'TAT')

        '''
        
        res = map(''.join, zip(*[iter(self._seq)] * 3))
        return tuple(map(lambda c: self.__class__(c), res)) if keep_class else tuple(res)
    

    def frames(self, complete=False):
        '''
        There are 3 reading frames on each strand of nucleotide (+1, +2, +3),
        and two strands (plus, minus) leading to complete 6 reading frames.
        The frames() method returns a dictionary of 3 reading frames by default, 
        or all 6 frames of the sequence if desired by user.

        >>> dna = DNA("AATTGGCCGCTTAAT")
        >>> dna.frames()
        {1: DNA(AATTGGCCGCTTAAT), 2: DNA(ATTGGCCGCTTAAT), 3: DNA(TTGGCCGCTTAAT)}
        >>> dna.frames(complete=True)
        {1: DNA(AATTGGCCGCTTAAT), 2: DNA(ATTGGCCGCTTAAT), 3: DNA(TTGGCCGCTTAAT), 
            -1: DNA(TTAACCGGCGAATTA), -2: DNA(TAACCGGCGAATTA), -3: DNA(AACCGGCGAATTA)}
        
        '''

        result = {}
        indices = (1, 2, 3, -1, -2, -3) if complete else (1, 2, 3)
        for frame in indices:
                if frame < 0:
                    result[frame] = self.c()[((frame * -1) - 1):]
                else:
                    result[frame] = self[(frame - 1):]
        return result


class DNA(AbstractNucleotide):
    '''
    DNA is the intended concrete class for all analysis using DNA sequences,
    which includes genome arithmetics and other DNA-related queries.
    '''
    
    def __init__(self, seq, valid_checker=Seq_data.valid_dna):      
        super().__init__(seq, valid_checker)


    def transcribe(self, mapper=Seq_data.default_transcribe):
        return super()._map(mapper, RNA)


    def t(self):
        return self.transcribe()


    def complement(self, mapper=Seq_data.default_complement_dna):
        return super().complement(mapper)


    def c(self, mapper=Seq_data.default_complement_dna):
        return self.complement(mapper)


    def reverse_complement(self, mapper=Seq_data.default_complement_dna):
        return self.complement(mapper)[::-1]


    def rc(self, mapper=Seq_data.default_complement_dna):
        return self.reverse_complement(mapper)


class RNA(AbstractNucleotide):
    def __init__(self, seq, valid_checker=Seq_data.valid_rna):
        super().__init__(seq, valid_checker)


    def reverse_transcribe(self, mapper=Seq_data.default_reverse_transcribe):
        return super()._map(mapper, DNA)
    

    def rt(self, mapper=Seq_data.default_reverse_transcribe):
        return self.reverse_transcribe(mapper)
    

    def complement(self, mapper=Seq_data.default_complement_rna):
        return super().complement(mapper)


    def c(self, mapper=Seq_data.default_complement_rna):
        return self.complement(mapper)


    def reverse_complement(self, mapper=Seq_data.default_complement_rna):
        return self.complement(mapper)[::-1]


    def rc(self, mapper=Seq_data.default_complement_rna):
        return self.reverse_complement(mapper)