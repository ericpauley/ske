package main

import (
	"bufio"
	"encoding/binary"
	"flag"
	"fmt"
	"io/ioutil"
	"os"
	"runtime"
	"sort"
	"strconv"
	"sync"
	"time"
	"unsafe"
)

func check(e error) {
	if e != nil {
		panic(e)
	}
}

const sectorexp = 14

var maxmem uint = 2048 * 1024 * 1024
var maxdisk uint = 10 * 1024 * 1024 * 1024
var maxCores uint
var minAbundance = 3
var counts countList

func calcSectors(f *os.File) []sector {
	kmers := make([]int, 0, 8*1024*1024)
	println("Calculating sectors")
	checked := 0
	scanned := scan(f, func(kmer kmer) bool {
		checked++
		if checked > 8*1024*1024 {
			return false
		}
		kmers = append(kmers, int(kmer.getPrefix()))
		return true
	}, false)
	sort.Ints(kmers)
	maxsize := uint(float64(maxmem/uint(unsafe.Sizeof(kmer{}))) * scanned)
	var sectors []sector
	var size uint
	var lastKmer int
	for _, kmer := range kmers {
		size++
		if size > maxsize {
			if lastKmer != kmer {
				s := sector{start: uint(lastKmer), end: uint(kmer)}
				sectors = append(sectors, s)
			}
			lastKmer = kmer
			size = 0
		}
	}
	sectors = append(sectors, sector{start: uint(lastKmer), end: uint(maxKmer.getPrefix() + 1)})
	return sectors
}

func saveChunks(oname string, counts countList, sectors *[]sector, start time.Time, sorted chan chan kmerlist, sjoin *sync.WaitGroup) {
	var outfiles []*outfile
	var writerJoin sync.WaitGroup
	for _, count := range counts {
		writerJoin.Add(1)
		println("Created count file ", count)
		name := oname + "." + strconv.Itoa(int(count))
		f, err := os.Create(name)
		check(err)
		ofile := outfile{file: f, stream: f, len: count}
		outfiles = append(outfiles, &ofile)
		ofile.channel = make(chan kmercount, 10000)

		go func() {
			for kc := range ofile.channel {
				fmt.Println("le", kc)
				fmt.Fprintln(ofile.stream, kc)
				/*kc.kmer.write(ofile.stream)
				b := make([]byte, 2)
				binary.LittleEndian.PutUint16(b, kc.count)
				ofile.stream.Write(b)*/
			}
			writerJoin.Done()
		}()
	}
	i := 0
	println("Ready to save sectors")
	for towrite := range sorted {
		var dcount uint32 = 1
		var kmers kmerlist
		if towrite != nil {
			println("Awaiting sector ", i)
			kmers = <-towrite
			println("Sector arrived ", len(kmers))
			i++
		} else {
			dcount = 0
			for i := 0; i < len(counts); i++ {
				kmers = append(kmers, maxKmer)
			}
			println("Recieved nil sector")
		}
		for j, kmer := range kmers {
			count := dcount
			len := kmer.autoRsh()
			for _, outfile := range outfiles {
				if len > outfile.len {
					kmer.rsh(uint(2 * (len - outfile.len)))
					len = outfile.len
				}
				if len == outfile.len {
					if outfile.count == 0 {
						outfile.current = kmer
					}
					if kmer == outfile.current && count != 0 {
						outfile.count += count
						break
					} else {
						fmt.Println("Overwrite", j, kmer)
						outfile.current, kmer = kmer, outfile.current
						outfile.count, count = count, outfile.count
						println(count)
						if int(count) >= minAbundance {
							v := kmercount{kmer: kmer, count: count}
							outfile.channel <- v
						}
					}
				}
			}
		}
	}
	for _, ofile := range outfiles {
		close(ofile.channel)
	}
	writerJoin.Wait()
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
		kmer.write(bufd)
	}
	bufd.Flush()
	pjoin.Done()
}

func processChunks(c chan kmerlist, sorted chan chan kmerlist, ojoin *sync.WaitGroup) {
	println("Process started")
	for kmers := range c {
		println("Sorting chunk!")
		output := make(chan kmerlist)
		sorted <- output
		sort.Sort(kmers)
		output <- kmers
	}
	ojoin.Done()
}

func revcmp() {
	f, err := os.Open("ecoli.count.31")
	reader := bufio.NewReaderSize(f, 4*1024*1024)
	check(err)
	b := make([]byte, 2)
	var counts kmercountlist
	for true {
		k := readKmer(reader).minimize(31)
		l, err := reader.Read(b)
		if err != nil {
			break
		}
		if l < 2 {
			break
		}
		count := binary.LittleEndian.Uint32(b)
		counts = append(counts, kmercount{k, count})
	}
	sort.Sort(counts)
	for _, count := range counts {
		fmt.Println(count.kmer, count.count)
	}
}

func main() {
	var foutput string
	flag.StringVar(&foutput, "out", "", "The output filename")
	flag.Var(&counts, "counts", "A comma separated list of kmer lengths to calculate")
	flag.UintVar(&maxmem, "maxmem", 2048, "Amount of memory allowed (MB)")
	flag.UintVar(&maxdisk, "maxdisk", 10, "Amount of disk usage allowed (GB)")
	flag.UintVar(&maxCores, "cores", uint(runtime.NumCPU()), "Number of CPU cores to use")
	flag.IntVar(&minAbundance, "min-abundance", 3, "Min number of occurences to be solid")
	flag.Parse()
	maxmem = maxmem * 1024 * 1024 / (maxCores + 2)
	maxdisk = maxdisk * 1024 * 1024 * 1024
	var finput = flag.Arg(0)
	if finput == "" {
		fmt.Println("Error: Must define an input file!")
		return
	}
	start := time.Now()
	oname := "ecoli.count"
	f, err := os.Open(finput)
	check(err)
	sectors := calcSectors(f)
	fmt.Printf("Created %d sectors\n", len(sectors))
	// toSort := make(chan sector)
	maxDiskSectors := maxdisk / maxmem
	if maxDiskSectors > 20 {
		maxDiskSectors = 20
	}
	passes := uint(len(sectors)-1)/maxDiskSectors + 1
	toSort := make(chan kmerlist)
	sorted := make(chan chan kmerlist, 100)
	var ojoin sync.WaitGroup
	var sjoin sync.WaitGroup
	for i := 0; uint(i) < maxCores; i++ {
		ojoin.Add(1)
		go processChunks(toSort, sorted, &ojoin)
	}
	sjoin.Add(1)
	go saveChunks(oname, counts, &sectors, start, sorted, &sjoin)
	for pass := uint(0); pass < passes; pass++ {
		fmt.Print("\rPerforming sectoring pass ", pass+1, " of ", passes)
		min := pass * maxDiskSectors
		max := (pass + 1) * maxDiskSectors
		if max > uint(len(sectors)) {
			max = uint(len(sectors))
		}
		toProcess := sectors[min:max]
		var pjoin sync.WaitGroup
		pjoin.Add(len(toProcess))
		var sectorMap [1 << 12]sector
		for i := range toProcess {
			toProcess[i].c = make(chan kmer, 100)
			go writeChunk(&toProcess[i], &pjoin)
			for j := toProcess[i].start; j < toProcess[i].end; j++ {
				sectorMap[j] = toProcess[i]
			}
		}
		lowest := maxKmer
		scan(f, func(kmer kmer) bool {
			if kmer.cmp(lowest) == -1 {
				lowest = kmer
			}
			if sectorMap[kmer.getPrefix()].c != nil {
				sectorMap[kmer.getPrefix()].c <- kmer
			}
			return true
		}, true)
		for _, s := range toProcess {
			close(s.c)
		}
		pjoin.Wait()
		for i := range toProcess {
			println("Reading sector")
			data := make(kmerlist, sectors[i].len)
			sectors[i].tmpfile.Seek(0, 0)
			reader := bufio.NewReaderSize(sectors[i].tmpfile, 4096*1024*1024)
			for j := 0; j < sectors[i].len; j++ {
				data[j] = readKmer(reader)
			}
			sectors[i].tmpfile.Close()
			err = os.Remove(sectors[i].tmpfile.Name())
			check(err)
			toSort <- data
		}
	}
	close(toSort)
	println("Waiting sort completion")
	ojoin.Wait()
	sorted <- nil
	println("waiting save completion")
	close(sorted)
	sjoin.Wait()
	fmt.Println("Counting took", time.Now().Sub(start), oname)
}
