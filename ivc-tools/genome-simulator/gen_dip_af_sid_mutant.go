/*----------------------------------------------------------------------------------
 Generate diploid mutant genomes based on SNPs and INDELs from VCF file.
 By Nam S. Vo, March 2015
 Generate diploid mutations (homozygous and heterozygous subs, indels) based on allele frequency.
 Usage: go run gen_dip_af_sid_mutant.go ref_file variant_file result_dir
 		E.g.: go run gen_dip_af_sid_mutant.go /data/nsvo/test-data/GRCh37_chr1/refs/GRCh37_chr1.fasta /data/nsvo/test-data/GRCh37_chr1/refs/TRIMMED.ALL.chr1.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf /data/nsvo/test-data/GRCh37_chr1/refs/dip_af_sid_mutant
 Output: - A simulated mutant genome:
		+ mutant_genome.fasta: mutant genome generated based on SNPs and INDELs from VCF file. Ignore other types of SVs when generating mutant genome
		+ (sub/ins/del)_same_ref.txt: set of variants of the mutant which are same as ref
		+ (sub/ins/del)_diff_ref.txt: set of variants of the mutant which are different from ref
		+ all_var.txt: set of all variants of the mutant
----------------------------------------------------------------------------------*/

package main

import (
	"bufio"
	"bytes"
	"fmt"
	"math/rand"
	"os"
	"path"
	"strconv"
	"time"
)

var IUPAC_nt_code = map[string]string{
	"AG": "R", "GA": "R",
	"CT": "Y", "TC": "Y",
	"GC": "S", "CG": "S",
	"AT": "W", "TA": "W",
	"GT": "K", "TG": "K",
	"AC": "M", "CA": "M",
}

func main() {
	if len(os.Args) < 4 {
		fmt.Println("Usage: go run gen_af_sid_mutant.go <input.fasta> <input.vcf> <result_dir>")
		return
	}
	fmt.Println("Reading reference...")
	chr_pos, chr_name, _ := GetGenome(os.Args[1])
	fmt.Println("Reading variant profile and generating mutations...")
	GetVarProf(os.Args[2], os.Args[3], chr_pos, chr_name)
	fmt.Println("Done!")
}

//--------------------------------------------------------------------------------------------------
// GetGenome gets reference genome from FASTA files.
//--------------------------------------------------------------------------------------------------
func GetGenome(file_name string) (chr_pos []int, chr_name [][]byte, seq []byte) {
	f, e := os.Open(file_name)
	if e != nil {
		panic(e)
	}
	defer f.Close()
	scanner := bufio.NewScanner(f)
	seq = make([]byte, 0)
	chr_pos = make([]int, 0)
	chr_name = make([][]byte, 0)
	var line []byte
	var sub_line [][]byte
	for scanner.Scan() {
		line = scanner.Bytes()
		if len(line) == 0 {
			continue
		}
		if line[0] == '>' {
			sub_line = bytes.Split(line, []byte(" "))
			chr_pos = append(chr_pos, len(seq))
			contig_name := make([]byte, len(sub_line[0][1:]))
			copy(contig_name, sub_line[0][1:])
			chr_name = append(chr_name, contig_name)
		} else {
			seq = append(seq, line...)
		}
	}
	return chr_pos, chr_name, seq
}

//--------------------------------------------------------------------------------------------------
// GetVarProf gets variant profile from VCF files.
//--------------------------------------------------------------------------------------------------
func GetVarProf(file_name, dir_name string, chr_pos []int, chr_name [][]byte) {

	var line, sline, sample, tmp_name, sub_info []byte
	var sub_line, sub_info_part [][]byte
	var i, j, pos, prev_pos int
	var sample_num, freq, hom_ref_freq, hom_var_freq, het_01_freq, het_10_freq float64
	var var_freq []float64
	var contig string

	rand.Seed(time.Now().UnixNano())

	if _, e := os.Stat(dir_name); os.IsNotExist(e) {
		os.Mkdir(dir_name, 0777)
	}
	mut_var_file, e := os.Create(path.Join(dir_name, path.Base(file_name)+".mutations.txt"))
	if e != nil {
		return
	}
	defer mut_var_file.Close()

	f, e := os.Open(file_name)
	if e != nil {
		panic(e)
	}
	defer f.Close()
	r := bufio.NewReader(f)
outer_loop:
	for {
		line, e = r.ReadBytes('\n')
		if e != nil {
			break
		}
		sline = bytes.Trim(line, "\n\r")
		if sline[0] == '#' || len(sline) == 0 {
			continue
		} else {
			sub_line = bytes.Split(sline, []byte("\t"))
			pos, _ = strconv.Atoi(string(sub_line[1]))
			if len(sub_line[3]) > 1 && prev_pos+len(sub_line[3]) > pos {
				continue
			}
			prev_pos = pos
			for _, sub_info = range bytes.Split(sub_line[7], []byte(";")) {
				sub_info_part = bytes.Split(sub_info, []byte("="))
				if bytes.Equal(sub_info_part[0], []byte("VT")) {
					if string(sub_info_part[1]) != "SNP" && string(sub_info_part[1]) != "INDEL" { //ignore other types of variants
						continue outer_loop
					}
				}
			}
			hom_ref, het_01, het_10, hom_var := 0, 0, 0, 0
			for _, sample = range sub_line[9:] {
				if sample[0] == '0' && sample[2] == '0' {
					hom_ref += 1
				}
				if sample[0] == '0' && sample[2] == '1' {
					het_01 += 1
				}
				if sample[0] == '1' && sample[2] == '0' {
					het_10 += 1
				}
				if sample[0] == '1' && sample[2] == '1' {
					hom_var += 1
				}
			}
			sample_num = float64(hom_ref + het_01 + het_10 + hom_var)
			hom_ref_freq = float64(hom_ref) / sample_num
			hom_var_freq = float64(hom_var) / sample_num
			het_01_freq = float64(het_01) / sample_num
			het_10_freq = float64(het_10) / sample_num
			var_freq = []float64{hom_ref_freq, hom_ref_freq + hom_var_freq, hom_ref_freq + hom_var_freq + het_01_freq, hom_ref_freq + hom_var_freq + het_01_freq + het_10_freq}
			//Get mutations
			freq = rand.Float64()
			//contig = string(sub_line[0])
			contig = "1"
			for i = 0; i < 4; i++ {
				if freq < var_freq[i] {
					if len(sub_line[3]) > 1 { //deletions
						for j = 1; j < len(sub_line[3]); j++ {
							if i == 3 { //het_10
								mut_var_file.WriteString(contig + "\t" + strconv.Itoa(pos+j) + "\t" + string(sub_line[3][j]) + "\t-\t1\n")
							} else if i == 2 { //het_01
								mut_var_file.WriteString(contig + "\t" + strconv.Itoa(pos+j) + "\t" + string(sub_line[3][j]) + "\t-\t2\n")
							} else if i == 1 { //hom_var
								mut_var_file.WriteString(contig + "\t" + strconv.Itoa(pos+j) + "\t" + string(sub_line[3][j]) + "\t-\t3\n")
							}
						}
					} else if len(sub_line[4]) > 1 { //insertions
						if i == 3 { //het_10
							mut_var_file.WriteString(contig + "\t" + string(sub_line[1]) + "\t-\t" + string(sub_line[4][1:]) + "\t1\n")
						} else if i == 2 { //het_01
							mut_var_file.WriteString(contig + "\t" + string(sub_line[1]) + "\t-\t" + string(sub_line[4][1:]) + "\t2\n")
						} else if i == 1 { //hom_var
							mut_var_file.WriteString(contig + "\t" + string(sub_line[1]) + "\t-\t" + string(sub_line[4][1:]) + "\t3\n")
						}
					} else if len(sub_line[3]) == 1 && len(sub_line[4]) == 1 { //substitutions
						if i == 3 { //het_10
							mut_var_file.WriteString(contig + "\t" + string(sub_line[1]) + "\t" + string(sub_line[3]) + "\t" + IUPAC_nt_code[string(sub_line[3])+string(sub_line[4])] + "\t1\n")
						} else if i == 2 { //het_01
							mut_var_file.WriteString(contig + "\t" + string(sub_line[1]) + "\t" + string(sub_line[3]) + "\t" + IUPAC_nt_code[string(sub_line[3])+string(sub_line[4])] + "\t2\n")
						} else if i == 1 { //hom_var
							mut_var_file.WriteString(contig + "\t" + string(sub_line[1]) + "\t" + string(sub_line[3]) + "\t" + string(sub_line[4]) + "\t3\n")
						}
					}
					break
				}
			}
			pos = -1
			for i, tmp_name = range chr_name {
				if string(sub_line[0]) == string(tmp_name) {
					pos = chr_pos[i]
					break
				}
			}
			if pos == -1 {
				fmt.Println("Contig", contig, "is missing from the reference")
			}
		}
	}
}
