package main

import (
	"bufio"
	"os"
)

func main() {
	input := bufio.NewScanner(os.Stdin)
	output := bufio.NewWriter(os.Stdout)
	push := false
	for input.Scan() {
		line := input.Bytes()
		if line[0] == '>' || line[0] == '@' {
			line[0] = '>'
			output.Write(line)
			output.WriteRune('\n')
			push = true
		} else if line[0] == '+' {
			push = false
		} else if push {
			output.Write(line)
			output.WriteRune('\n')
		}
	}
}
