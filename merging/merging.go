package main

import (
	"flag"
	"fmt"
	"runtime/debug"
	"strconv"
	"time"

	"github.com/ericpauley/dna"
	"github.com/ericpauley/go-hdf5"
)

const readLimit = 10000

var minAbundance uint = 3
var minsize uint = 8
var maxsize uint = 30

func check(e error) {
	if e != nil {
		panic(e)
	}
}

func kmerReader(table *hdf5.Table, reqs chan tableRead, size int, records int) chan []dna.Kmer {
	c := make(chan []dna.Kmer, 10)
	go func() {
		for j := 0; j < records; j += readLimit {
			toFetch := readLimit
			if j+toFetch > records {
				toFetch = records - j
			}
			if toFetch == 0 {
				return
			}
			req := tableRead{table, toFetch, make(chan []dna.Kmer)}
			reqs <- req
			result := <-req.future
			if len(result) > 0 && result[0].Length == 0 {
				for i := range result {
					result[i].Normalize(uint32(size))
				}
			}
			tosend := result[:0]
			for i := range result {
				if result[i].Length >= uint32(minsize) {
					result[i].Truncate(uint32(maxsize))
					tosend = append(tosend, result[i])
				}
			}
			c <- tosend
		}
		close(c)
	}()
	return c
}

func mergeStreams(first chan []dna.Kmer, second chan []dna.Kmer) chan []dna.Kmer {
	c := make(chan []dna.Kmer, 1)
	buf := make([]dna.Kmer, 0, readLimit)
	var current dna.Kmer
	go func() {
		var fVals, sVals []dna.Kmer
		frun, srun, ok := true, true, true
		for frun || srun {
			var toPush dna.Kmer
			if frun && len(fVals) == 0 {
				fVals, ok = <-first
				if !ok {
					frun = false
				}
			}
			if srun && len(sVals) == 0 {
				sVals, ok = <-second
				if !ok {
					srun = false
				}
			}
			if frun && srun {
				if fVals[0].Cmp(sVals[0]) < 0 {
					toPush = fVals[0]
					fVals = fVals[1:]
				} else {
					toPush = sVals[0]
					sVals = sVals[1:]
				}
			} else if frun && !srun {
				toPush = fVals[0]
				fVals = fVals[1:]
			} else if srun && !frun {
				toPush = sVals[0]
				sVals = sVals[1:]
			} else {
				break
			}
			if current.Cmp(toPush) != 0 {
				if current.Count > 0 {
					buf = append(buf, current)
				}
				current = toPush
			} else {
				current.Count += toPush.Count
			}
			if len(buf) >= readLimit {
				c <- buf
				buf = make([]dna.Kmer, 0, readLimit)
			}
		}
		c <- buf
		close(c)
	}()
	return c
}

type tableRead struct {
	table  *hdf5.Table
	num    int
	future chan []dna.Kmer
}

type tableWrite struct {
	table  *hdf5.Table
	data   []minimerCount
	future chan error
}

func streamKmers(reads chan tableRead, writes chan tableWrite, done chan bool) {
	for {
		select {
		case read := <-reads:
			result := make([]dna.Kmer, read.num, readLimit)
			if len(result) > 0 {
				read.table.Next(&result)
			}
			read.future <- result
		case write := <-writes:
			select {
			case write.future <- write.table.Append(&write.data):
			default:
			}
		case <-done:
			done <- true
			return
		}

	}
}

type minimerCount struct {
	mmer  dna.Minimer
	count uint32
}

type output struct {
	table   *hdf5.Table
	file    *hdf5.File
	current minimerCount
	buffer  []minimerCount
}

func main() {
	flag.UintVar(&minAbundance, "min-abundance", 1, "Min number of occurences to be solid")
	flag.UintVar(&minsize, "min-size", 8, "Min kmer size to count")
	flag.UintVar(&maxsize, "max-size", 30, "Max kmer size to count")
	flag.Parse()
	names := flag.Args()
	if len(names) == 0 {
		fmt.Println("Error: Must define an input file!")
		return
	}
	var kmersources []chan []dna.Kmer
	reads := make(chan tableRead)
	writes := make(chan tableWrite, 2)
	streamWait := make(chan bool)
	go streamKmers(reads, writes, streamWait)
	for _, name := range names {
		h5, err := hdf5.OpenFile(name, hdf5.F_ACC_RDONLY)
		check(err)
		group, err := h5.OpenGroup("dsk/solid")
		var size int
		if err != nil {
			group, err = h5.OpenGroup("partials")
			check(err)
		} else {
			size = 31
		}
		num, err := group.NumObjects()
		check(err)
		for i := uint(0); i < num; i++ {
			name, err := group.ObjectNameByIndex(i)
			check(err)
			dset, _ := group.OpenTable(name)
			records, _ := dset.NumPackets()
			kmersources = append(kmersources, kmerReader(dset, reads, size, records))
		}
	}
	for len(kmersources) > 1 {
		kmersources = append(kmersources[2:], mergeStreams(kmersources[0], kmersources[1]))
	}
	outputs := make([]output, maxsize+1)
	for i := maxsize; i >= minsize; i-- {
		h5, err := hdf5.CreateFile("merged"+strconv.Itoa(int(i))+".h5", hdf5.F_ACC_TRUNC)
		check(err)
		table, err := h5.CreateTableFrom("kmers", minimerCount{}, 1<<20, -1)
		check(err)
		outputs[i].table = table
		outputs[i].file = h5
		outputs[i].buffer = make([]minimerCount, 0, readLimit)
	}
	i := 0
	for _, kmersource := range kmersources {
		for countlist := range kmersource {
			for _, kmer := range countlist {
				i++
				if i%1000000 == 0 {
					debug.FreeOSMemory()
					fmt.Println(i, kmer.ToRaw(), kmer.Length, kmer.Count)
				}
				for l := kmer.Length; l >= uint32(minsize); l-- {
					kmer.Truncate(l)
					mmer := kmer.ToRaw()
					if mmer.Cmp(outputs[l].current.mmer) == 0 {
						outputs[l].current.count += kmer.Count
					} else {
						if outputs[l].current.count >= uint32(minAbundance) {
							outputs[l].buffer = append(outputs[l].buffer, outputs[l].current)
						}
						outputs[l].current = minimerCount{mmer, kmer.Count}
						if len(outputs[l].buffer) >= readLimit {
							writes <- tableWrite{outputs[l].table, outputs[l].buffer, make(chan error)}
							outputs[l].buffer = make([]minimerCount, 0, readLimit)
						}
					}
				}
			}
		}
	}
	close(reads)
	close(writes)
	streamWait <- true
	<-streamWait
	time.Sleep(time.Second)
	for i := maxsize; i >= minsize; i-- {
		outputs[i].table.Close()
		outputs[i].file.Flush(hdf5.F_SCOPE_GLOBAL)
		outputs[i].file.Close()
	}

	println(i)
}
