package main

import (
	"bufio"
	"bytes"
	"fmt"
	"os"
	"sync"
	"time"

	"github.com/ericpauley/dna"
)

func parser(ch chan []byte, join *sync.WaitGroup, tocall dna.Kmerhandler, running *bool, min int, max int) {
	defer join.Done()
	for s := range ch {
		var kmer dna.Kmer
		var pushed int
		for i := len(s) - 1; i >= 0; i-- {
			if s[i]>>3&1 == 0 {
				bp := uint64(s[i] >> 1 & 3)
				kmer.Push(bp)
				pushed++
				if pushed == max {
					if !tocall(kmer) {
						*running = false
					}
				}
			} else {
				pushed = 0
				kmer = dna.Kmer{}
			}
		}
		if pushed < max && pushed >= min {
			if !tocall(kmer) {
				*running = false
			}
		}
	}
}

func scan(f *os.File, tocall dna.Kmerhandler, min int, max int, verbose bool) float64 {
	f.Seek(0, 0)
	start := time.Now()
	fi, err := f.Stat()
	check(err)
	size := fi.Size()
	scanner := bufio.NewScanner(bufio.NewReader(f))
	running := true

	c := make(chan []byte)
	var pjoin sync.WaitGroup
	pjoin.Add(1)
	for i := 0; i < 1; i++ {
		go parser(c, &pjoin, tocall, &running, min, max)
	}
	index := 0
	var current []byte
	line := 0
	for scanner.Scan() {
		line++
		b := scanner.Bytes()
		if bytes.IndexByte(b, '>') != -1 {
			index++
			if index%10000 == 0 {
				pos, err := f.Seek(0, 1)
				check(err)
				if verbose {
					fmt.Print("\r", scanner.Text(), len(current), pos*100/size, time.Since(start), time.Since(start).Nanoseconds()/1000000*size/pos, "          ")
				}
			}
			c <- current
			current = make([]byte, 0)
		} else {
			current = append(current, b...)
		}
		if !running {
			break
		}
	}
	c <- current
	close(c)
	pjoin.Wait()
	pos, err := f.Seek(0, 1)
	check(err)
	return float64(pos) / float64(size)
}
