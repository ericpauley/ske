package main

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

type sector struct {
	start int
	end   int
}

func swap(x *uint64, y *uint64) {
	temp := *y
	*y = *x
	*x = temp
}

func partition(array []uint64, p uint, q uint, pivotLocation uint) uint {

	pivot := array[pivotLocation]
	swap(&array[pivotLocation], &array[q])
	i := p
	for j := p; j < q; j++ {
		if array[j] <= pivot {
			swap(&array[i], &array[j])
			i++
		}
	}
	swap(&array[q], &array[i])
	return i
}

func quicksort(array []uint64, start uint, end uint) {
	if start < end {
		pivot := (end + start) / 2
		r := partition(array, start, end, pivot)
		if r > start {
			quicksort(array, start, r-1)
		}
		quicksort(array, r+1, end)
	}
}

func check(e error) {
	if e != nil {
		panic(e)
	}
}

const sectorexp = 14
const maxmem = 512 * 1024 * 1024

func parser(ch chan string, join *sync.WaitGroup, tocall kmerhandler, running *bool) {
	defer join.Done()
	for s := range ch {
		var kmer kmer
		kmer.init()
		var pushed bool
		for len, c := range s {
			switch {
			case c == 'A' || c == 'a':
				pushed = true
			case c == 'C' || c == 'c':
				kmer.push(1)
				pushed = true
			case c == 'T' || c == 't':
				kmer.push(2)
				pushed = true
			case c == 'G' || c == 'g':
				kmer.push(3)
				pushed = true
			default:
				kmer.init()
				pushed = false
			}
			if pushed {
				if !tocall(kmer) {
					*running = false
				}
			}
		}
	}
}

func scan(f *os.File, tocall kmerhandler, verbose bool) float64 {
	f.Seek(0, 0)
	start := time.Now()
	fi, err := f.Stat()
	check(err)
	size := fi.Size()
	scanner := bufio.NewScanner(bufio.NewReaderSize(f, 1024*1024))

	running := true

	c := make(chan string, 2000)
	var pjoin sync.WaitGroup
	pjoin.Add(10)
	for i := 0; i < 10; i++ {
		go parser(c, &pjoin, tocall, &running)
	}
	index := 0
	current := ""

	for scanner.Scan() {
		if strings.HasPrefix(scanner.Text(), ">") {
			index++
			if index%1000000 == 0 {
				pos, err := f.Seek(0, 1)
				check(err)
				if verbose {
					fmt.Println(scanner.Text(), len(current), pos*100/size, time.Since(start), time.Since(start).Nanoseconds()/1000000*size/pos)
				}
			}
			c <- current
			current = ""
		} else {
			current += scanner.Text()
		}
		if !running {
			break
		}
	}
	close(c)
	pjoin.Wait()
	pos, err := f.Seek(0, 1)
	check(err)
	fmt.Println("Scan took ", time.Since(start))
	return float64(pos) / float64(size)
}

func calcSectors(f *os.File) ([]sector, []int) {
	kmers := make(kmerlist, 0, 16*1024*1024)
	println("Calculating sectors")
	checked := 0
	scanned := scan(f, func(kmer kmer) bool {
		checked++
		if checked > 16*1024*1024 {
			return false
		}
		kmers = append(kmers, kmer)
		return true
	}, false)
	sort.Sort(kmers)

	for i := range sectors {
		sectors[i] = int(float64(sectors[i]) / scanned * 1.1)
	}
	size := 0
	sum := 0
	slices := 1
	for i := range sectors {
		if size+sectors[i] > maxmem/8 {
			slices++
			size = 0
		}
		size += sectors[i]
		sum += sectors[i]
	}
	lastIndex := 0
	size = 0
	dsectors := make([]sector, slices)
	sectornum := 0
	for i := range sectors {
		if size+sectors[i] > maxmem/8 || size > sum/slices {
			dsectors[sectornum] = sector{start: lastIndex, end: i}
			println("Sector", sectornum, "from", lastIndex, "to", i, "size", size*8)
			lastIndex = i
			size = 0
			sectornum++
		}
		size += sectors[i]
	}
	dsectors[sectornum] = sector{start: lastIndex, end: 1 << sectorexp}
	println("Sector", sectornum, "from", lastIndex, "to", 1<<sectorexp, "size", size*8)
	return dsectors, sectors
}

type kmercount struct {
	kmer  uint64
	count uint16
}

type outfile struct {
	file    *os.File
	stream  io.Writer
	len     int
	current uint64
	count   uint16
	channel chan kmercount
}

func saveChunks(oname string, counts []int, sorted chan []uint64, start time.Time) {
	var outfiles []*outfile
	for _, count := range counts {
		name := oname + "." + strconv.Itoa(count)
		f, err := os.Create(name)
		check(err)
		ofile := outfile{file: f, stream: bufio.NewWriterSize(f, 1024*1024), len: count}
		outfiles = append(outfiles, &ofile)
		ofile.channel = make(chan kmercount, 10000)
		go func() {
			for kc := range ofile.channel {
				binary.Write(ofile.stream, binary.LittleEndian, kc.kmer)
				binary.Write(ofile.stream, binary.LittleEndian, kc.count)
			}
		}()
	}
	i := 0
	for c := range sorted {
		if i%100 == 0 {
			println("Outputting chunk", i, i*100/(1<<sectorexp), time.Since(start).Seconds()*float64(1<<sectorexp)/float64(i))
		}
		i++
		for _, kmer := range c {
			len := 31
			var count uint16 = 1
			for kmer&3 == 0 {
				kmer >>= 2
				len--
			}
			kmer >>= 2
			for _, outfile := range outfiles {
				if len > outfile.len {
					kmer >>= uint(2 * (len - outfile.len))
					len = outfile.len
				}
				if len == outfile.len {
					if kmer == outfile.current || outfile.count == 0 {
						outfile.count++
						break
					} else {
						outfile.current, kmer = kmer, outfile.current
						outfile.count, count = count, outfile.count
						outfile.channel <- kmercount{kmer: kmer, count: count}
					}
				}
			}
		}
	}
}

func processChunks(c chan []uint64, sorted chan []uint64) {
	i := 0
	for sector := range c {
		i++
		sort.Sort(kmerlist(sector))
		sorted <- sector
	}
}

func main() {
	start := time.Now()
	fname := "ecoli.fasta"
	oname := "ecoli.count"
	counts := []int{31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1}
	f, err := os.Open(fname)
	check(err)
	sectors, sizes := calcSectors(f)
	fmt.Printf("Created %d sectors\n", len(sectors))
	c := make(chan []uint64, 1000)
	sorted := make(chan []uint64, 1000)
	go processChunks(c, sorted)
	go saveChunks(oname, counts, sorted, start)
	for i, sector := range sectors {
		fmt.Printf("Creating sector #%d\n", i)
		chunks := make([][]uint64, sector.end-sector.start)
		for i := range chunks {
			chunks[i] = make([]uint64, 0, sizes[i+sector.start])
		}
		scan(f, func(kmer uint64) bool {
			cnum := int(kmer >> (64 - sectorexp))
			if sector.start <= cnum && cnum < sector.end {
				chunks[cnum-sector.start] = append(chunks[cnum-sector.start], kmer)
			}
			return true
		}, true)
		for i, chunk := range chunks {
			if len(chunk) > sizes[i+sector.start] {
				println("chunk went over", i+sector.start)
			}
			c <- chunk
		}
	}
}
