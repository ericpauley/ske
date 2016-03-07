package main

import (
	"bufio"
	"fmt"
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
			fmt.Fprintln(output, line)
			push = true
		} else if line[0] == '+' {
			push = false
		} else if push {
			fmt.Fprintln(output, line)
		}
	}
}
