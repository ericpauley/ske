package main

import (
	"encoding/binary"
	"fmt"
	"io"
	"math/rand"
	"os"
	"strconv"
	"strings"
)

const kmerwords uint = 1
const maxCount = 32*kmerwords - 1

type kmer [kmerwords]uint64

func (kmer *kmer) push(bp uint64) {
	var carry = bp << 62
	pushed := (kmer[kmerwords-1]&3 == 3)
	for i := range kmer {
		carry, kmer[i] = kmer[i]<<62, kmer[i]>>2|carry
	}
	if pushed {
		kmer[kmerwords-1] |= 3
	}
}

func (kmer *kmer) rsh(count uint) {
	var carry uint64
	for i := range kmer {
		carry, kmer[i] = kmer[i]<<(64-count), kmer[i]>>count|carry
	}
}

func (kmer *kmer) lsh(count uint) {
	var carry uint64
	for i := int(kmerwords - 1); i >= 0; i-- {
		carry, kmer[i] = kmer[i]>>(64-count), kmer[i]<<count|carry
	}
}

func (kmer *kmer) autoRsh() uint {
	len := maxCount
	for kmer[kmerwords-1]&3 != 3 {
		kmer.rsh(2)
		len--
	}
	kmer.rsh(2)
	return len
}

func (kmer *kmer) init() {
	for i := range kmer {
		kmer[i] = 0
	}
	kmer[0] = uint64(3) << 62
}

func (kmer *kmer) cmp(rhs kmer) int {
	for i := range kmer {
		if kmer[i] < rhs[i] {
			return -1
		} else if kmer[i] > rhs[i] {
			return 1
		}
	}
	return 0
}

func (kmer *kmer) getPrefix() uint {
	return uint(kmer[0] >> (64 - 12))
}

func (kmer *kmer) write(w io.Writer) {
	var b [kmerwords * 8]byte
	for i := 0; uint(i) < kmerwords; i++ {
		binary.LittleEndian.PutUint64(b[i*8:(i+1)*8], uint64(kmer[i]))
	}
	w.Write(b[:])
}

func (k kmer) minimize(len int) kmer {
	var output kmer
	orig := k
	for i := 0; i < len; i++ {
		bp := (orig[kmerwords-1] + 2) & 3
		orig.rsh(2)
		output.lsh(2)
		output[kmerwords-1] |= bp
	}
	switch k.cmp(output) {
	case -1:
		return k
	case 1:
		return output
	default:
		return k
	}
}

func readKmer(r io.Reader) kmer {
	var b = make([]byte, kmerwords*8)
	r.Read(b)
	var kmer kmer
	for i := 0; uint(i) < kmerwords; i++ {
		kmer[i] = binary.LittleEndian.Uint64(b[i*8 : (i+1)*8])
	}
	return kmer
}

var maxKmer kmer

type kmerhandler func(kmer) bool

type kmerlist []kmer

func (a kmerlist) Len() int           { return len(a) }
func (a kmerlist) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a kmerlist) Less(i, j int) bool { return a[i].cmp(a[j]) == -1 }

func qsort(a kmerlist) kmerlist {
	if len(a) < 2 {
		return a
	}

	left, right := 0, len(a)-1

	// Pick a pivot
	pivotIndex := rand.Int() % len(a)

	// Move the pivot to the right
	a[pivotIndex], a[right] = a[right], a[pivotIndex]

	// Pile elements smaller than the pivot on the left
	for i := range a {
		if a[i].cmp(a[right]) == -1 {
			a[i], a[left] = a[left], a[i]
			left++
		}
	}

	// Place the pivot after the last smaller element
	a[left], a[right] = a[right], a[left]

	// Go down the rabbit hole
	qsort(a[:left])
	qsort(a[left+1:])

	return a
}

type countList []uint

func (i *countList) String() string {
	return fmt.Sprint(*i)
}

func (i *countList) Set(value string) error {
	// If we wanted to allow the flag to be set multiple times,
	// accumulating values, we would delete this if statement.
	// That would permit usages such as
	//	-deltaT 10s -deltaT 15s
	// and other combinations.
	for _, dt := range strings.Split(value, ",") {
		count, err := strconv.Atoi(dt)
		if err != nil {
			return err
		}
		*i = append(*i, uint(count))
	}
	return nil
}

func (i *countList) Max() uint {
	var max uint
	for _, j := range *i {
		if j > max {
			max = j
		}
	}
	return max
}

type sector struct {
	start   uint
	end     uint
	tmpfile *os.File
	c       chan kmer
	len     int
}

type kmercount struct {
	kmer  kmer
	count uint16
}

type kmercountlist []kmercount

func (a kmercountlist) Len() int           { return len(a) }
func (a kmercountlist) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a kmercountlist) Less(i, j int) bool { return a[i].kmer.cmp(a[j].kmer) == -1 }

type outfile struct {
	file    *os.File
	stream  io.Writer
	len     uint
	current kmer
	count   uint16
	channel chan kmercount
}

type empty struct{}

type semaphore chan empty

func (semaphore semaphore) acquire() {
	semaphore <- empty{}
}

func (semaphore semaphore) release() {
	<-semaphore
}

func init() {
	for i := range maxKmer {
		maxKmer[i] = ^uint64(0)
	}
}
