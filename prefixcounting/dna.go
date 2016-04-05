package main

import (
	"bufio"
	"flag"
	"fmt"
	"io/ioutil"
	"os"
	"runtime"
	"runtime/debug"
	"strconv"
	"strings"
	"sync"
	"time"
	"unsafe"

	"github.com/ericpauley/dna"
	"github.com/sbinet/go-hdf5"
)

func check(e error) {
	if e != nil {
		panic(e)
	}
}

const sorters = 1

var maxmem uint = 2048 * 1024 * 1024
var maxCores uint
var minAbundance = 1
var counts countList
var minsize = 8
var maxsize = 30
var maxrecords uint

type sector struct {
	tmpfile *os.File
	c       chan dna.Kmer
	len     int
	sorted  chan dna.Mmerlist
	toSort  chan dna.Mmerlist
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

func calcSectors(f *os.File) ([]*sector, int) {
	println("Calculating sectors")
	checked := 0
	scanned := scan(f, func(kmer dna.Kmer) bool {
		checked++
		if checked > 1024 {
			return false
		}
		return true
	}, minsize, maxsize, false)
	total := int(float64(maxsize-minsize+1) * float64(checked) / scanned)
	sectors := 4
	println(total / sectors)
	sectorbits := 2
	for uint(total/sectors) > maxrecords {
		sectors <<= 1
		sectorbits++
	}
	sectorSlice := make([]*sector, sectors)
	for i := 0; i < sectors; i++ {
		sectorSlice[i] = &sector{nil, make(chan dna.Kmer, 100), 0, make(chan dna.Mmerlist), make(chan dna.Mmerlist)}
	}
	return sectorSlice, sectorbits
}

func saveChunks(oname string, counts countList, sectors []*sector, sjoin *sync.WaitGroup) {
	h5, err := hdf5.CreateFile(oname, hdf5.F_ACC_TRUNC)
	check(err)
	group, _ := h5.CreateGroup("partials")
	println("Ready to save sectors")
	start := time.Now()
	for snum, sector := range sectors {
		kmers := <-sector.sorted
		fmt.Println("Received sector", len(kmers), time.Now().Sub(start))
		ostart := time.Now()
		table, err := group.CreateTableFrom(strconv.Itoa(snum), dna.Kmer{}, 1<<20, -1)
		check(err)
		todump := make(dna.Kmerlist, 0, 10000)
		var current dna.Kmer
		for _, mmer := range kmers {
			kmer := mmer.ToKmer()
			if current.Cmp(kmer) == 0 {
				current.Count++
			} else {
				if current.Length >= uint32(minAbundance) {
					todump = append(todump, current)
				}
				current = kmer
			}
			if len(todump) >= 10000 {
				table.Append(&todump)
				todump = todump[:0]
			}
		}
		kmers = nil
		debug.FreeOSMemory()
		table.Append(&todump)
		table.Close()
		h5.Flush(hdf5.F_SCOPE_GLOBAL)
		fmt.Println("Saved", time.Now().Sub(ostart))
		start = time.Now()
	}
	h5.Flush(hdf5.F_SCOPE_GLOBAL)
	h5.Close()
	sjoin.Done()
}

func writeChunk(s *sector, pjoin *sync.WaitGroup) {
	t, err := ioutil.TempFile("", "kmer")
	check(err)
	s.tmpfile = t
	bufd := bufio.NewWriterSize(t, 4*1024*1024)
	check(err)
	for kmer := range s.c {
		s.len++
		kmer.Write(bufd)
	}
	bufd.Flush()
	pjoin.Done()
}

func processChunks(c chan *sector, ojoin *sync.WaitGroup) {
	println("Process started")
	start := time.Now()
	for sector := range c {
		kmers := <-sector.toSort
		fmt.Println("Sorting chunk!", time.Now().Sub(start))
		kmers.Sort(0)
		sector.sorted <- kmers
		kmers = nil
		debug.FreeOSMemory()
		start = time.Now()
	}
	ojoin.Done()
}

func main() {
	fmt.Println(os.TempDir())
	var foutput string
	flag.StringVar(&foutput, "out", "", "The output filename")
	flag.UintVar(&maxmem, "maxmem", 2048, "Amount of memory allowed (MB)")
	flag.UintVar(&maxCores, "cores", uint(runtime.NumCPU()), "Number of CPU cores to use")
	flag.IntVar(&minAbundance, "min-abundance", 1, "Min number of occurences to be solid")
	flag.IntVar(&minsize, "min-size", 8, "Min kmer size to count")
	flag.IntVar(&maxsize, "max-size", 30, "Max kmer size to count")
	flag.Parse()
	maxmem = maxmem * 1024 * 1024
	maxrecords = maxmem / (sorters + 2) / uint(unsafe.Sizeof(dna.Minimer{}))
	var finput = flag.Arg(0)
	if finput == "" {
		fmt.Println("Error: Must define an input file!")
		return
	}
	start := time.Now()
	if foutput == "" {
		parts := strings.Split(finput, ".")
		if len(parts) > 1 {
			parts = parts[:len(parts)-1]
		}
		foutput = strings.Join(parts, ".") + ".partials.h5"
	}
	f, err := os.Open(finput)
	check(err)
	sectors, sectorbits := calcSectors(f)
	fmt.Printf("Created %d sectors\n", len(sectors))
	var pjoin sync.WaitGroup
	pjoin.Add(len(sectors))
	for i := range sectors {
		go writeChunk(sectors[i], &pjoin)
	}
	scan(f, func(kmer dna.Kmer) bool {
		hash := kmer.ToRaw()[len(kmer.Kmer)-1]
		for i := 0; i < minsize*2-sectorbits; i++ {
			hash ^= hash >> 1
		}
		sectors[hash&uint64(len(sectors)-1)].c <- kmer
		return true
	}, minsize, maxsize, true)
	for _, s := range sectors {
		close(s.c)
	}
	pjoin.Wait()
	toSort := make(chan *sector)
	var ojoin sync.WaitGroup
	var sjoin sync.WaitGroup
	for i := 0; i < sorters; i++ {
		ojoin.Add(1)
		go processChunks(toSort, &ojoin)
	}
	sjoin.Add(1)
	go saveChunks(foutput, counts, sectors, &sjoin)
	for _, sector := range sectors {
		println("Reading sector")
		data := make(dna.Mmerlist, 0, sector.len*(maxsize-minsize+1))
		println(cap(data))
		sector.tmpfile.Seek(0, 0)
		reader := bufio.NewReaderSize(sector.tmpfile, 4*1024*1024)
		for j := 0; j < sector.len; j++ {
			k := dna.ReadKmer(reader)
			for k.Length >= uint32(minsize) {
				data = append(data, k.ToMini())
				k.Cut()
			}
		}
		println(cap(data), len(data))
		sector.tmpfile.Close()
		err = os.Remove(sector.tmpfile.Name())
		check(err)
		toSort <- sector
		sector.toSort <- data
	}
	close(toSort)
	println("Waiting sort completion")
	ojoin.Wait()
	println("waiting save completion")
	sjoin.Wait()
	fmt.Println("Counting took", time.Now().Sub(start), foutput)
}
