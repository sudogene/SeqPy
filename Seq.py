import Seq_data
import warnings


class Seq(str):
    '''
    Seq (Sequence) is an abstract subclass of string with added
    functionality for biological sequences. The Seq class is not
    intended to be instantiated, and it is recommended for users
    to work with concrete classes like DNA, RNA, Protein.

    >>> string = "ATGC"
    >>> Seq(string)
    Seq(ATGC)
    '''

    def __init__(self, seq):
        self = seq.upper()


    def __repr__(self):
        LIMIT = 21
        if len(self) > LIMIT:
            return f'{self.__class__.__name__}({self[:LIMIT]}...{self[-3:]})'
        return f'{self.__class__.__name__}({self})'


    def __eq__(self, other):
        if self.__class__ == other.__class__:
            return str(self) == str(other)
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


    def _map(self, mapper, seqclass=None, window=1):
        '''
        The generic map function for all sequence mapping, allows for a 
        robust abstract method for all variations of mapping (mapper), 
        and the return sequence subclass (seqclass).
        '''
        if not seqclass:
            seqclass = self.__class__
        try:
            if window > 1:
                groups = [self[i:i+window] for i in range(0, len(self), window)]
                filtered = list(filter(lambda ele: len(ele) == window, groups))
                if len(groups) != len(filtered):
                    warnings.formatwarning = lambda msg, *args, **kwargs: f'{msg}\n'
                    warnings.warn(f'Some elements do not fit window={window} and are ignored.')
                mapped = [mapper[group] for group in filtered]
            else:
                mapped = [mapper[ele] for ele in self]
            new_seq = ''.join(mapped)
            return seqclass(new_seq)
        except KeyError:
            raise KeyError('The provided mapper does not contain some bases in the sequence.')
    

    def reverse(self):
        return self.__class__(self[::-1])
    

    def count(self, pattern):
        if not self._class_check(pattern):
            raise TypeError('Invalid input type')
        return self.count(str(pattern))


    def count_overlap(self, pattern):
        if not self._class_check(pattern):
            raise TypeError('Invalid input type')
        
        pattern = str(pattern)
        count = start = 0
        while True:
            start = self.find(pattern, start) + 1
            if start > 0:
                count += 1
            else:
                return count


class AbstractNucleotide(Seq):
    '''
    AbstractNucleotide is the abstract class for nucleotide sequences 
    including DNA and RNA. It has parent functions for finding complements, 
    generating reading frames and codons.
    '''

    def __init__(self, seq, valid_checker):
        '''
        Each nucleotide type has a set of valid bases in the 
        valid_checker argument. Most subclasses of AbstractNucleotide 
        have their own default sets;
        e.g. {A, T, G, C, N} for DNA and {A, U, G, C, N} for RNA.

        Users can also insert their own set if desired.
        '''

        for base in seq:
            if base not in valid_checker:
                raise ValueError('Invalid sequence')       
        super().__init__(seq)


    def complement(self, mapper):
        '''
        Returns complementary sequence
        '''

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
    

    def kmers(self, k):
        '''
        Returns a list of all k-mers, subsequence of nucleotide of length k.

        >>> DNA('CAATCCAAC').get_kmer(5)
        ['CAATC', 'AATCC', 'ATCCA', 'TCCAA', 'CCAAC']

        '''
        
        res = [self[i: i+k] for i in range(len(self)-k+1)]
        return res


    def codons(self):
        '''
        A codon is a nucleotide with 3 bases, read during translation event.
        Returns a tuple of codon strings by default, else a tuple of codon nucleotides.

        >>> DNA("AGGCCTGGGGTTTAATTAGCGTAGCGTTTACGATAT").codons()
        ['AGG', 'CCT', 'GGG', 'GTT', 'TAA', 'TTA', 'GCG', 'TAG', 'CGT', 'TTA', 'CGA', 'TAT']

        '''
        
        res = map(''.join, zip(*[iter(self)] * 3))
        return list(res)
    

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
    

    def translate(self):
        return self.transcribe().translate()


class RNA(AbstractNucleotide):
    ''' 
    RNA is the intended concrete class for all analysis using RNA sequences.
    '''

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
    

    def translate(self, mapper=Seq_data.default_translate):
        return super()._map(mapper, Protein, 3)


class Protein(Seq):
    '''
    Protein sequence
    '''
    def __init__(self, seq, valid_checker=Seq_data.valid_protein):
        for base in seq:
            if base not in valid_checker:
                raise ValueError('Invalid sequence')       
        super().__init__(seq)
    

