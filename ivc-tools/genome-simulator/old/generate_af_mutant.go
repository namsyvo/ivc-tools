///////////////////////////////////////////////////////////////////////////////////////
// Mutate genome using SNPs and INDELs (vcf file) with a reference genome (fasta file)
// Fall 2013 - Quang
// Nam S. Vo, Sep 2014
// Modified by Nam S. Vo, Jan 2015: generate mutations based on allele frequency
///////////////////////////////////////////////////////////////////////////////////////

package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"bytes"
	"sort"
	"time"
	"math/rand"
	"path"
	"runtime"
	"log"
	"math"
)

type SNPProfile struct {
	Profile    [][]byte
	AlleleFreq []float32
}

var NEW_VAR_PROB = []float64{0.00}
var STD_BASES = []byte{'A', 'C', 'G', 'T'}

//Global variable for memory profiling
var Memstats = new(runtime.MemStats)

func PrintMemStats(mesg string) {
    runtime.ReadMemStats(Memstats)
    log.Printf(mesg + "\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f",
		Memstats.Alloc, Memstats.TotalAlloc, Memstats.Sys, Memstats.HeapAlloc, Memstats.HeapSys,
		float64(Memstats.Alloc)/(math.Pow(1024, 3)), float64(Memstats.Sys)/(math.Pow(1024, 3)))
}

func main(){
	if len(os.Args) < 4 {
		fmt.Println("Usage: go run genome-mutate.go <input.fasta> <input.vcf> <result_dir>")
		return
	}
	fmt.Println("Reading FASTA file...")
	header, seq := ReadFASTA(os.Args[1])
	fmt.Println("Reading VCF file...")
	SNP_Prof := ReadVCF(os.Args[2])
	for _, prob := range NEW_VAR_PROB {
		fmt.Println("Generating mutated genome " + strconv.FormatFloat(prob, 'f', 5, 32))
		mutate, known_variants, new_snps, new_indels := GenerateMutate(seq, SNP_Prof, prob)
		fmt.Println("Saving everything...")
		SaveMutations(path.Join(os.Args[3], "new-mutate-" + strconv.FormatFloat(prob, 'f', 5, 32)), header, mutate, known_variants, new_snps, new_indels)		
	}
}

//--------------------------------------------------------------------------------------------------
// Generate mutated genome
//--------------------------------------------------------------------------------------------------
func GenerateMutate(seq []byte, SNP_Prof map[int]SNPProfile, new_var_prob float64) ([]byte, map[int][]byte, map[int]byte, map[int][]byte) {
	mutate := make([]byte, 0)
	known_variants := make(map[int][]byte) //known snps and indels
	new_snps := make(map[int]byte) // new snps (at known locations)
	new_indels := make(map[int][]byte) // new indels (at known locations)

	snp_pos := make([]int, 0, len(SNP_Prof))
	for key, _ := range SNP_Prof {
		snp_pos = append(snp_pos, key)

	}
	sort.Ints(snp_pos)
	fmt.Println("# SNPs", len(snp_pos))

	var pos, i, idx, snp_idx, change_idx int
	var b byte

	PrintMemStats("Memstats after generating new SNPs")

	mutate = append(mutate, seq[ : snp_pos[0]]...)

	rand.Seed(time.Now().UnixNano()) // takes the current time in nanoseconds as the seed
	snp_idx = snp_pos[0]
	idx = rand.Intn(len(SNP_Prof[snp_idx].Profile))
	known_variants[snp_idx] = SNP_Prof[snp_idx].Profile[idx]
	if string(known_variants[snp_idx]) != "." || string(known_variants[snp_idx]) != "<DEL>" {
		mutate = append(mutate, SNP_Prof[snp_idx].Profile[idx]...)
	}
	var snp_prof [][]byte
	var is_indel bool
	var allele_freq []float32
	sub_ins_count, del_count := 0, 0
	for snp_idx = 1; snp_idx < len(snp_pos); snp_idx++ {
		pos = snp_pos[snp_idx]
		snp_prof = SNP_Prof[pos].Profile
		allele_freq = SNP_Prof[pos].AlleleFreq
		mutate = append(mutate, seq[snp_pos[snp_idx - 1] + 1 : pos]...)
		if rand.Float32() < float32(new_var_prob) {
			var snps []byte
			var indels [][]byte
			is_indel = false
			for idx = 0; idx < len(snp_prof); idx++ {
				if len(snp_prof[idx]) == 1 && string(snp_prof[idx][0]) != "."{
					snps = append(snps, snp_prof[idx][0])
				} else {
					indels = append(indels, snp_prof[idx])
					is_indel = true
				}
			}
			if !is_indel {
				for {
					i = rand.Intn(4)
					b = STD_BASES[i]
					if bytes.IndexByte(snps, b) == -1 {
						mutate = append(mutate, b)
						new_snps[pos] = b
						break
					}
				}
			} else {
				if string(indels[0]) != "." {
					indel := make([]byte, len(indels[0]))
					copy(indel, indels[0])
					change_idx = 1 + rand.Intn(len(indel) - 1) //do not change the first base
					for {
						i = rand.Intn(4)
						b = STD_BASES[i]
						if bytes.IndexByte(indel[change_idx : change_idx + 1], b) == -1 {
							indel[change_idx] = b
							mutate = append(mutate, indel...)
							new_indels[pos] = indel
							break
						}
					}
				} else {
					new_indels[pos] = indels[0]
				}
			}
		} else {
			if string(snp_prof[1]) != "." {
				if rand.Float32() < allele_freq[1] {
					known_variants[pos] = snp_prof[1]
					mutate = append(mutate, snp_prof[1]...)
				} else {
					known_variants[pos] = snp_prof[0]
					mutate = append(mutate, snp_prof[0]...)
				}
				sub_ins_count++
			} else {
				known_variants[pos] = snp_prof[1]
				del_count++
			}
		}
	}
	fmt.Println("SUB_INS, DEL count", sub_ins_count, del_count)
	mutate = append(mutate, seq[snp_pos[len(snp_pos) - 1] + 1 : ]...)
	PrintMemStats("Memstats after generating new alleles")

	return mutate, known_variants, new_snps, new_indels
}

//--------------------------------------------------------------------------------------------------
// Save mutated genome, mutation profiles to file
//--------------------------------------------------------------------------------------------------
func SaveMutations(dir_name string, header, mutate []byte, known_variants map[int][]byte, new_snps map[int]byte, new_indels map[int][]byte) {
	if _, err := os.Stat(dir_name); os.IsNotExist(err) {
	    os.Mkdir(dir_name, 0777)
	}

	ref_file, err := os.Create(path.Join(dir_name, "mutant_genome.fasta"))
    if err != nil {
        // handle the error here
        return
    }
    defer ref_file.Close()
	ref_file.Write(header)
	ref_file.WriteString("\n")
	ref_file.Write(mutate)

	kv_file, err := os.Create(path.Join(dir_name, "known_variants.txt"))
	if  err != nil {
        return
    }
	defer kv_file.Close()

	var known_variant_pos []int
	for key, _ := range known_variants {
		 known_variant_pos = append(known_variant_pos, key)
	}
	sort.Ints(known_variant_pos)

	for _, pos := range known_variant_pos {
		snp_pos := strconv.Itoa(pos)
		_, err := kv_file.WriteString(snp_pos + "\t" + string(known_variants[pos]) + "\n");
		if err != nil {
            fmt.Println(err)
            break
        }
	}

	ns_file, err := os.Create(path.Join(dir_name, "new_snps.txt"))
	if  err != nil {
        return
    }
	defer ns_file.Close()

	var new_snp_pos []int
	for key, _ := range new_snps {
		 new_snp_pos = append(new_snp_pos, key)
	}
	sort.Ints(new_snp_pos)

	for _, pos := range new_snp_pos {
		snp_pos := strconv.Itoa(pos)
		_, err := ns_file.WriteString(snp_pos + "\t" + string(new_snps[pos]) + "\n");
		if err != nil {
            fmt.Println(err)
            break
        }
	}

	nid_file, err := os.Create(path.Join(dir_name, "new_indels.txt"))
	if  err != nil {
        return
    }
	defer nid_file.Close()

	var new_indel_pos []int
	for key, _ := range new_indels {
		 new_indel_pos = append(new_indel_pos, key)
	}
	sort.Ints(new_indel_pos)

	for _, pos := range new_indel_pos {
		snp_pos := strconv.Itoa(pos)
		_, err := nid_file.WriteString(snp_pos + "\t" + string(new_indels[pos]) + "\n");
		if err != nil {
            fmt.Println(err)
            break
        }
	}
	PrintMemStats("Memstats after writing mutate")
}

//--------------------------------------------------------------------------------------------------
// Read FASTA files
//--------------------------------------------------------------------------------------------------
func ReadFASTA(sequence_file string) ([]byte, []byte) {
	f, err := os.Open(sequence_file)
	if err != nil {
		fmt.Printf("%v\n", err)
		os.Exit(1)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	scanner.Scan()
	line := scanner.Bytes()
	var header = make([]byte, len(line))
	copy(header, line)
	var seq []byte
	for scanner.Scan() {
		line = scanner.Bytes()
		seq = append(seq, line...)
	}
	return header, seq
}

//--------------------------------------------------------------------------------------------------
// Read VCF files
//--------------------------------------------------------------------------------------------------
func ReadVCF(sequence_file string) map[int]SNPProfile {
	SNP_Prof := make(map[int]SNPProfile)
	f, err := os.Open(sequence_file)
	if err != nil {
		fmt.Printf("%v\n", err)
		os.Exit(1)
	}
	defer f.Close()

	var line, sub_info []byte
	var alt, info, sub_info_part [][]byte
	var pos int
	var p float32

	rand.Seed(time.Now().UnixNano()) // takes the current time in nanoseconds as the seed
	data := bufio.NewReader(f)
	for {
		line, err = data.ReadBytes('\n')
		if err != nil {
			fmt.Printf("%v\n",err)
			break
		}
		if bytes.Equal(line[0:1], []byte("#")) {
			//Do something here later
		} else {
			tmp := SNPProfile{}
			sub_line, _ := SplitN(line, []byte("\t"), 9)
			pos, _ = strconv.Atoi(string(sub_line[1]))

			//Take REF bases if they are deletions
			//if len(sub_line[3]) > 1 {
			//	tmp.Profile = append(tmp.Profile, sub_line[3])
			//}

			//Take only location which has Sub or Ins
			//Put into SNP Prof only ALT bases, not REF bases
			if len(sub_line[3]) == 1 {
				if len(alt) > 1 {
					fmt.Println("Not a biallelic mutation", pos, len(alt))
				}
				tmp.Profile = append(tmp.Profile, sub_line[3])
				p = rand.Float32()
				tmp.AlleleFreq = append(tmp.AlleleFreq, p)
				alt = bytes.Split(sub_line[4], []byte(","))
				info = bytes.Split(sub_line[7], []byte(";"))
				for _, sub_info = range info {
					sub_info_part = bytes.Split(sub_info, []byte("="))
					if bytes.Equal(sub_info_part[0], []byte("AF")) {
						tmp_p, _ := strconv.ParseFloat(string(sub_info_part[1]), 32)
						p = float32(tmp_p)
						break
					}
				}
				for i := 0; i < len(alt); i++ {
					if bytes.Equal(alt[i], []byte("<DEL>")) {
						tmp.Profile = append(tmp.Profile, []byte("."))
					} else {
						tmp.Profile = append(tmp.Profile, alt[i])
					}
					tmp.AlleleFreq = append(tmp.AlleleFreq, p)
				}
				SNP_Prof[pos - 1] = tmp // append SNP at pos
			}
		}
	}
	PrintMemStats("Memstats after reading SNP profile")
	return SNP_Prof
}

//--------------------------------------------------------------------------------------------------
//Memory-efficient string split function
//--------------------------------------------------------------------------------------------------
func SplitN(s, sep []byte, n int) ([][]byte, int) {
	first_idx, sep_idx := 0, 0
	sep_num := 0
	t := make([][]byte, 0)
	for first_idx < len(s) {
		sep_idx = bytes.Index(s[first_idx : ], sep)
		if sep_idx != -1 {
			sep_num++
			tmp := make([]byte, first_idx + sep_idx - first_idx)
			copy(tmp, s[first_idx : first_idx + sep_idx])
			t = append(t, tmp)
			if sep_num == n {
				return t, sep_num
			}
			first_idx = first_idx + sep_idx + 1
		} else {
			return t, sep_num
		}
	}
	return t, sep_num
}