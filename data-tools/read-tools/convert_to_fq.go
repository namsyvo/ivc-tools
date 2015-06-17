/*

Convert the simulated reads to .fq format

EX:
GGAGAGATTTGCAAAAAGATAGGATGAGGTTGCATGAATTATCACCTTTTACAAAAGGGCAGTCAATCTGTGTTTTCCAAAAAAAAAAATCAAATCAAAA 1 13578847 1 50

@r1
GGAGAGATTTGCAAAAAGATAGGATGAGGTTGCATGAATTATCACCTTTTACAAAAGGGCAGTCAATCTGTGTTTTCCAAAAAAAAAAATCAAATCAAAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

*/

package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"
//	"strconv"
	"bytes"
//	"unicode/utf8"
)

//import "time"
//import "log"

func usage() {
	fmt.Println("Usage: go run convert_to_fq.go <input.txt>")
	fmt.Println("Example: go run convert_to_fq.go test.txt > test.fq")
}

func main() {

	if len(os.Args) != 2 {
		usage()
		os.Exit(0)
	}

	//Q = (ERR_RATE == 0.0)? 'I' : (int)(-10.0 * log(ERR_RATE) / log(10.0) + 0.499) + 33;
	n := 0
	f,err := os.Open(os.Args[1])
    if err != nil{
        fmt.Printf("%v\n",err)
        os.Exit(1)
    }

    defer f.Close()
    br := bufio.NewReader(f)

    for {
	    // line , isPrefix, err := br.ReadLine()
	    line , err := br.ReadString('\n')
	    if err != nil {
			//fmt.Printf("%v\n",err)
			break
		} else {
			n++			
			sline := string(line)
			split := strings.Split(sline, " ");
			fmt.Println(fmt.Sprintf("@r%d", n))
			fmt.Printf("%s\n",split[0])			
			fmt.Println("+")
			var buffer bytes.Buffer

		    for i := 0; i < len(split[0]); i++ {
		        buffer.WriteString("I")
		    }

    		fmt.Println(buffer.String())			
		}		    
	}	
}