/*

Extract the alternative in a .vcf file from GATK

*/

package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"
//	"strconv"
//	"bytes"
//	"sort"
)

//import "time"
//import "math/rand"

func usage() {
	fmt.Println("Usage: go run buildindex.go <input.vcf> <chr_name>")
	fmt.Println("Example: go run buildindex.go test.vcf soy_bean")
}

func main(){
	if len(os.Args) != 3 {
		usage()
		os.Exit(0)
	}

	vcfRead(os.Args[1])	
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
		line , err := br.ReadString('\n')
		if err != nil {
			//fmt.Printf("%v\n",err)
			break
		}
		if line[0]==byte('#') {
			fmt.Printf("%s",line)
		} else {
			sline := string(line)
			split := strings.Split(sline, "\t")
			fmt.Print(strings.Replace(sline, split[0], os.Args[2], 1))
		}
	}
    return array
}