/*

Extract the alternative in a .vcf file from GATK

*/

package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"
	"strconv"
//	"bytes"
//	"sort"
)

//import "time"
//import "math/rand"

func usage() {
	fmt.Println("Usage: go run buildindex.go <input.vcf> ")
	fmt.Println("Example: go run buildindex.go test.vcf")
}

func main(){
	if len(os.Args) != 2 {
		usage()
		os.Exit(0)
	}

	SNP_profile := vcfRead(os.Args[1])
	
	for key, v := range SNP_profile {
		fmt.Println(key, " ", v)
	}
	
}

func vcfRead(sequence_file string) map[int]string {
	array := make(map[int]string)
	f,err := os.Open(sequence_file)
    if err != nil{
        fmt.Printf("%v\n",err)
        os.Exit(1)
    }

    defer f.Close()
    br := bufio.NewReader(f)
    //byte_array := bytes.Buffer{}

	for{
		line , isPrefix, err := br.ReadLine()
		if err != nil || isPrefix{
			//fmt.Printf("%v\n",err)
			break
		}		
		if line[0]==byte('#') {
			//fmt.Printf("%s \n",line)
		} else {
			sline := string(line)
			split := strings.Split(sline, "\t");
			//fmt.Printf("%s %s %s\n", split[1], split[3], split[4])
			pos, _ := strconv.ParseInt(split[1], 10, 64)
			array[int(pos)]= split[4] // append SNP at pos
		}
	}
    return array
}