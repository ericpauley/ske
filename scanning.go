package main

import (
	"bufio"
	"fmt"
	"os"
	"sync"
	"time"
)

func parser(ch chan []byte, join *sync.WaitGroup, tocall kmerhandler, running *bool) {
	defer join.Done()
	for s := range ch {
		var kmer kmer
		kmer.init()
		var pushed bool
		for _, c := range s {
			kmer.push(uint64(c >> 1 & 3))
			if kmer[0] == 3 {
				fmt.Println(s)
			}
			pushed = true
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

	c := make(chan []byte)
	var pjoin sync.WaitGroup
	pjoin.Add(1)
	for i := 0; i < 1; i++ {
		go parser(c, &pjoin, tocall, &running)
	}
	index := 0
	var current []byte

	for scanner.Scan() {
		b := scanner.Bytes()
		if b[0] == '>' {
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
			if len(current) == 0 {
				current = b
			} else {
				current = append(current, b...)
			}
		}
		if !running {
			break
		}
	}
	close(c)
	pjoin.Wait()
	pos, err := f.Seek(0, 1)
	check(err)
	return float64(pos) / float64(size)
}
