package dna

import (
	"encoding/binary"
	"io"
	"math/rand"
)

const kmerwords uint = 1
const maxCount = 32*kmerwords - 1

type Minimer [kmerwords]uint64

/*Kmer represents a kmer and its count */
type Kmer struct {
	Kmer   Minimer `value`
	Count  uint32  `abundance`
	Length uint32  `length`
}

// Push a new base pair to the kmer
func (kmer *Kmer) Push(bp uint64) {
	var carry = bp << 62
	for i := range kmer.Kmer {
		carry, kmer.Kmer[i] = kmer.Kmer[i]<<62, kmer.Kmer[i]>>2|carry
	}
	kmer.Length++
}

func (kmer *Kmer) Cut() {
	kmer.Kmer.Lsh(2)
	kmer.Length--
}

// Rsh right shifts the kmer count bits
func (kmer *Minimer) Rsh(count uint32) {
	var carry uint64
	for i := range kmer {
		carry, kmer[i] = kmer[i]<<(64-count), kmer[i]>>count|carry
	}
}

// Lsh left shifts the kmer count bits
func (kmer *Minimer) Lsh(count uint32) {
	var carry uint64
	for i := len(kmer) - 1; i >= 0; i-- {
		carry, kmer[i] = kmer[i]>>(64-count), kmer[i]<<count|carry
	}
}

// Normalize takes a LSB aligned Kmer and converts it to a MSB aligned kmer
func (kmer *Kmer) Normalize(length uint32) {
	kmer.Kmer.Lsh(uint32(64*len(kmer.Kmer)) - 2*length)
	kmer.Length = length
}

// Truncate the kmer to length base pairs
func (kmer *Kmer) Truncate(length uint32) {
	if kmer.Length > length {
		for i := range kmer.Kmer {
			if uint32(i) > kmer.Length/32 {
				kmer.Kmer[i] = 0
			} else if uint32(i) == kmer.Length/32 {
				kmer.Kmer[i] = kmer.Kmer[i] & ^(^uint64(0) >> (2 * (length - uint32(32*i))))
			}
		}
		kmer.Length = length
	}
}

func (kmer *Kmer) ToRaw() Minimer {
	out := kmer.Kmer
	out.Rsh(uint32(64*kmerwords) - 2*kmer.Length)
	return out
}

// Cmp compares two Kmers
func (kmer Kmer) Cmp(rhs Kmer) int {
	mcmp := kmer.Kmer.Cmp(rhs.Kmer)
	if mcmp == 0 {
		return int(kmer.Length) - int(rhs.Length)
	}
	return mcmp
}

// Cmp compares two Kmers
func (kmer Minimer) Cmp(rhs Minimer) int {
	for i := range kmer {
		if kmer[i] < rhs[i] {
			return -1
		} else if kmer[i] > rhs[i] {
			return 1
		}
	}
	return 0
}

func (kmer Kmer) ToMini() Minimer {
	mmer := kmer.Kmer
	mmer[kmer.Length/64] |= (3 << (62 - (kmer.Length%32)*2))
	return mmer
}

func (kmer Minimer) ToKmer() Kmer {
	var out Kmer
	out.Length = uint32(kmerwords*32 - 1)
	for kmer[kmerwords-1]&3 == 0 {
		kmer.Rsh(2)
		out.Length--
	}
	kmer.Rsh(2)
	out.Kmer = kmer
	out.Normalize(out.Length)
	out.Count = 1
	return out
}

// GetPrefix returns the kmer prefix
func (kmer *Kmer) GetPrefix() uint {
	return uint(kmer.Kmer[0] >> (64 - 12))
}

func (kmer Minimer) Write(w io.Writer) {
	var b [kmerwords * 8]byte
	for i := 0; uint(i) < kmerwords; i++ {
		binary.LittleEndian.PutUint64(b[i*8:(i+1)*8], uint64(kmer[i]))
	}
	w.Write(b[:])
}

func (kmer Kmer) Write(w io.Writer) {
	kmer.ToMini().Write(w)
}

func ReadMinimer(r io.Reader) Minimer {
	var b = make([]byte, kmerwords*8)
	r.Read(b)
	var kmer Minimer
	for i := 0; uint(i) < kmerwords; i++ {
		kmer[i] = binary.LittleEndian.Uint64(b[i*8 : (i+1)*8])
	}
	return kmer
}

func ReadKmer(r io.Reader) Kmer {
	return ReadMinimer(r).ToKmer()
}

var maxKmer Kmer

type Kmerhandler func(Kmer) bool

type Kmerlist []Kmer

func (a Kmerlist) Len() int           { return len(a) }
func (a Kmerlist) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a Kmerlist) Less(i, j int) bool { return a[i].Cmp(a[j]) == -1 }

type Mmerlist []Minimer

func (a Mmerlist) Sort(depth int) {
	if len(a) < 2 {
		return
	}

	left, right := 0, len(a)-1

	// Pick a pivot
	pivotIndex := rand.Int() % len(a)

	// Move the pivot to the right
	a[pivotIndex], a[right] = a[right], a[pivotIndex]

	// Pile elements smaller than the pivot on the left
	for i := range a {
		if a[i].Cmp(a[right]) < 0 {
			a[i], a[left] = a[left], a[i]
			left++
		}
	}

	// Place the pivot after the last smaller element
	a[left], a[right] = a[right], a[left]

	// Go down the rabbit hole
	if depth < 5 {
		fin := make(chan bool)
		go func() {
			a[:left].Sort(depth + 1)
			fin <- true
		}()
		go func() {
			a[left+1:].Sort(depth + 1)
			fin <- true
		}()
		<-fin
		<-fin
	} else {
		a[:left].Sort(depth + 1)
		a[left+1:].Sort(depth + 1)
	}
}

func (a Mmerlist) Len() int           { return len(a) }
func (a Mmerlist) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a Mmerlist) Less(i, j int) bool { return a[i].Cmp(a[j]) == -1 }

/*type kmercountlist []Kmer

func (a kmercountlist) Len() int           { return len(a) }
func (a kmercountlist) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a kmercountlist) Less(i, j int) bool { return a[i].Cmp(a[j]) == -1 }

/*type outfile struct {
	file    *os.File
	stream  io.Writer
	len     uint
	current kmer
	count   uint32
	channel chan kmercount
}*/

type empty struct{}

type semaphore chan empty

func (semaphore semaphore) acquire() {
	semaphore <- empty{}
}

func (semaphore semaphore) release() {
	<-semaphore
}

var BlankMmer Minimer

func init() {
	for i := range maxKmer.Kmer {
		maxKmer.Kmer[i] = ^uint64(0)
		BlankMmer = Kmer{}.ToMini()
	}
}
