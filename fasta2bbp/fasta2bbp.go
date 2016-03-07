package main

import (
	"bufio"
	"encoding/binary"
	"os"
)

const (
	a = 0
	c = 1 << 62
	t = 2 << 62
	g = 3 << 62
)

func main() {
	var block uint64 = 1
	var pushed uint8
	dumpBuffer := make([]byte, 8)
	input := bufio.NewScanner(os.Stdin)
	output := bufio.NewWriter(os.Stdout)
	for input.Scan() {
		line := input.Bytes()
		if len(line) == 0 || line[0] == '>' || line[0] == '@' {
			if pushed != 0 {
				binary.LittleEndian.PutUint64(dumpBuffer, block)
				output.Write(dumpBuffer)
				pushed = 0
				block = c
			}
		} else if line[0] == '+' {
			if pushed != 0 {
				binary.LittleEndian.PutUint64(dumpBuffer, block)
				output.Write(dumpBuffer)
				pushed = 0
				block = c
			}
			input.Scan()
		} else {
			for _, b := range line {
				switch b {
				case 'A', 'a':
					block = block>>2 | a
				case 'C', 'c':
					block = block>>2 | c
				case 'T', 't':
					block = block>>2 | t
				case 'G', 'g':
					block = block>>2 | g
				default:
					binary.LittleEndian.PutUint64(dumpBuffer, block)
					output.Write(dumpBuffer)
					pushed = 0
					block = c
					continue
				}
				pushed++
				if pushed == 31 {
					binary.LittleEndian.PutUint64(dumpBuffer, block)
					output.Write(dumpBuffer)
					pushed = 0
					block = t
				}
			}
		}
	}
	output.Flush()
}
