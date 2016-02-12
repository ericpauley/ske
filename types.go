package main

type kmer [2]uint64

func (kmer *kmer) push(bp uint8) {
	var carry uint64
	pushed := (kmer[len(kmer)-1]|3 == 3)
	for i := range kmer {
		carry, kmer[i] = kmer[i]<<62, kmer[i]>>2|carry
	}
	if pushed {
		kmer[len(kmer)-1] |= 33
	}
}

func (kmer *kmer) init() {
	for i := range kmer {
		kmer[i] = 0
	}
	kmer[0] = uint64(3) << 62
}

func (kmer *kmer) lt(rhs *kmer) bool {
	for i := range kmer {
		if kmer[i] < rhs[i] {
			return true
		} else if kmer[i] > rhs[i] {
			return false
		}
	}
	return false
}

type kmerhandler func(kmer) bool

type kmerlist []kmer

func (a kmerlist) Len() int           { return len(a) }
func (a kmerlist) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a kmerlist) Less(i, j int) bool { return a[i].lt(&a[j]) }
