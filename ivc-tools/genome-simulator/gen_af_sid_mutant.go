/*----------------------------------------------------------------------------------
 Generate mutant genome based on SNPs and INDELs from VCF file
 By Nam S. Vo, March 2015
 Generate mutations (subs, indels) based on allele frequency
 Usage: go run generate_af_sid_mutant.go ref_file variant_file result_dir
 		E.g.: go run gen_af_sid_mutant.go /data/nsvo/test-data/GRCh37_chr1/refs/GRCh37_chr1.fasta /data/nsvo/test-data/GRCh37_chr1/refs/TRIMMED.ALL.chr1.phase1_release_v3.20101123.snps_indels_svs.genotypes.diffcontigname.vcf /data/nsvo/test-data/GRCh37_chr1/refs/af_sid_mutant
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
	"sort"
	"strconv"
	"time"
)

type VarProf struct {
	Variant [][]byte
	AleFreq []float32
	VarType []byte
}

func main() {
	if len(os.Args) < 4 {
		fmt.Println("Usage: go run gen_af_sid_mutant.go <input.fasta> <input.vcf> <result_dir>")
		return
	}
	fmt.Println("Reading reference...")
	chr_pos, chr_name, seq := GetGenome(os.Args[1])
	fmt.Println("Reading variant profile...")
	var_prof := GetVarProf(os.Args[2], chr_pos, chr_name)
	fmt.Println("Generating mutant genome and corresponding variant profile...")
	mutant, mut_var, all_mut_var := GenerateMutant(seq, var_prof)
	fmt.Println("Saving mutant gennome and its variant profile...")
	SaveMutant(os.Args[3], chr_pos, chr_name, mutant, mut_var, all_mut_var)
	fmt.Println("Done!")
}

//--------------------------------------------------------------------------------------------------
// Generate mutant genome
//--------------------------------------------------------------------------------------------------
func GenerateMutant(seq []byte, var_prof map[int]VarProf) ([]byte, []map[int][][]byte, map[int][][]byte) {

	mutant := make([]byte, 0)
	mut_var := make([]map[int][][]byte, 6)
	for i := 0; i < 6; i++ {
		mut_var[i] = make(map[int][][]byte)
	}
	all_mut_var := make(map[int][][]byte)

	snp_pos := make([]int, 0, len(var_prof))
	for key, _ := range var_prof {
		snp_pos = append(snp_pos, key)

	}
	sort.Ints(snp_pos)
	fmt.Println("Total number of variants:", len(snp_pos))

	var pos, idx int

	rand.Seed(time.Now().UnixNano())
	mutant = append(mutant, seq[:snp_pos[0]]...)

	var variant [][]byte
	var alefreq []float32
	var vartype []byte
	sv_num, ol_var_num := 0, 0
	for idx = 0; idx < len(snp_pos)-1; idx++ {
		pos = snp_pos[idx]
		variant = var_prof[pos].Variant
		alefreq = var_prof[pos].AleFreq
		vartype = var_prof[pos].VarType
		if len(vartype) != 0 && string(vartype) != "SNP" && string(vartype) != "INDEL" { //ignore other types of variants
			mutant = append(mutant, seq[snp_pos[idx]:snp_pos[idx+1]]...)
			sv_num++
		} else {
			if rand.Float32() < alefreq[1] { //Take ALT bases
				if len(variant[0]) == 1 { //SNPs or Insertions
					if len(variant[1]) == 1 { //SNPs
						mut_var[0][pos] = variant
					} else { //Insertions
						mut_var[1][pos] = variant
					}
					all_mut_var[pos] = variant
					mutant = append(mutant, variant[1]...)
					mutant = append(mutant, seq[snp_pos[idx]+1:snp_pos[idx+1]]...)
				} else if len(variant[0]) > 1 { //Deletions
					mut_var[2][pos] = variant
					all_mut_var[pos] = variant
					mutant = append(mutant, variant[1]...)
					//Ignore SNP profile that overlap with the selected deletion
					anchor_pos := snp_pos[idx] + len(variant[0])
					for snp_pos[idx+1] < anchor_pos {
						idx++
						ol_var_num++
					}
					mutant = append(mutant, seq[anchor_pos:snp_pos[idx+1]]...)
				}
			} else { //Take REF bases
				if len(variant[0]) == 1 { //SNPs or Insertions
					if len(variant[1]) == 1 { //SNPs
						mut_var[3][pos] = append(mut_var[3][pos], variant[0])
						mut_var[3][pos] = append(mut_var[3][pos], variant[0])
					} else { //Insertions
						mut_var[4][pos] = append(mut_var[4][pos], variant[0])
						mut_var[4][pos] = append(mut_var[4][pos], variant[0])
					}
					all_mut_var[pos] = append(all_mut_var[pos], variant[0])
					all_mut_var[pos] = append(all_mut_var[pos], variant[0])
					mutant = append(mutant, variant[0]...)
					mutant = append(mutant, seq[snp_pos[idx]+1:snp_pos[idx+1]]...)
				} else if len(variant[0]) > 1 { //Deletions
					mut_var[5][pos] = append(mut_var[5][pos], variant[0])
					mut_var[5][pos] = append(mut_var[5][pos], variant[0])
					all_mut_var[pos] = append(all_mut_var[pos], variant[0])
					all_mut_var[pos] = append(all_mut_var[pos], variant[0])
					mutant = append(mutant, variant[0]...)
					//Ignore SNP profile that overlap with the selected deletion
					anchor_pos := snp_pos[idx] + len(variant[0])
					for snp_pos[idx+1] < anchor_pos {
						idx++
						ol_var_num++
					}
					mutant = append(mutant, seq[anchor_pos:snp_pos[idx+1]]...)
				}
			}
		}
	}
	mutant = append(mutant, seq[snp_pos[idx]:]...) //ignore the last SNP at this moment
	fmt.Println("Total, sub_diff_ref, ins_diff_ref, del_diff_ref, sub_same_ref, ins_same_ref, del_same_ref, SV, OL_VAR:")
	fmt.Println(len(var_prof), len(mut_var[0]), len(mut_var[1]), len(mut_var[2]), len(mut_var[3]), len(mut_var[4]), len(mut_var[5]), sv_num, ol_var_num)

	return mutant, mut_var, all_mut_var
}

//--------------------------------------------------------------------------------------------------
// Save mutant genome, mutation profiles to file
//--------------------------------------------------------------------------------------------------
func SaveMutant(dir_name string, chr_pos []int, chr_name [][]byte, mutant []byte, mut_var []map[int][][]byte, all_mut_var map[int][][]byte) {
	if _, e := os.Stat(dir_name); os.IsNotExist(e) {
		os.Mkdir(dir_name, 0777)
	}

	f, e := os.Create(path.Join(dir_name, "mutant_genome.fasta"))
	if e != nil {
		// handle the error here
		return
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	for i, pos := range chr_pos {
		w.WriteString(">" + string(chr_name[i]) + "\t" + strconv.Itoa(pos) + "\n")
		if i < len(chr_pos)-1 {
			w.Write(mutant[pos:chr_pos[i+1]])
		} else {
			w.Write(mutant[pos:])
		}
		w.WriteString("\n")
	}
	w.Flush()

	file_names := []string{"sub_diff_ref.txt", "ins_diff_ref.txt", "del_diff_ref.txt", "sub_same_ref.txt", "ins_same_ref.txt", "del_same_ref.txt"}
	for i, file_name := range file_names {
		SaveVarProf(dir_name, file_name, mut_var[i])
	}
	SaveVarProf(dir_name, "mut_var.txt", all_mut_var)
}

func SaveVarProf(dir_name, file_name string, mut_var map[int][][]byte) {
	mut_var_file, e := os.Create(path.Join(dir_name, file_name))
	if e != nil {
		return
	}
	defer mut_var_file.Close()

	var mut_var_key []int
	for key, _ := range mut_var {
		mut_var_key = append(mut_var_key, key)
	}
	sort.Ints(mut_var_key)
	_, e = mut_var_file.WriteString("#POS\tREF\tALT\n")
	for _, pos := range mut_var_key {
		snp_pos := strconv.Itoa(pos)
		_, e = mut_var_file.WriteString(snp_pos + "\t" + string(mut_var[pos][0]) + "\t" + string(mut_var[pos][1]) + "\n")
		if e != nil {
			fmt.Println(e)
			break
		}
	}
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
func GetVarProf(file_name string, chr_pos []int, chr_name [][]byte) map[int]VarProf {
	var_prof := make(map[int]VarProf)
	f, e := os.Open(file_name)
	if e != nil {
		panic(e)
	}
	defer f.Close()

	var line, sline, info, sub_info, tmp_af, tmp_name []byte
	var sub_line, sub_info_part, info_arr [][]byte
	var i, pos, var_pos int
	var alt_prob float32
	var tmp_p float64
	var af []float32
	r := bufio.NewReader(f)
	for {
		line, e = r.ReadBytes('\n')
		if e != nil {
			break
		}
		sline = bytes.Trim(line, "\n\r")
		if sline[0] == '#' || len(sline) == 0 {
			continue
		} else {
			var_prof_elem := VarProf{}
			sub_line = bytes.SplitN(sline, []byte("\t"), 9)

			ref := make([]byte, len(sub_line[3]))
			copy(ref, sub_line[3])
			var_prof_elem.Variant = append(var_prof_elem.Variant, ref)

			alt := make([]byte, len(sub_line[4]))
			copy(alt, sub_line[4])
			alt_arr := bytes.Split(alt, []byte(","))
			info = make([]byte, len(sub_line[7]))
			copy(info, sub_line[7])
			info_arr = bytes.Split(sub_line[7], []byte(";"))
			var vt []byte
			af = make([]float32, 0)
			for _, sub_info = range info_arr {
				sub_info_part = bytes.Split(sub_info, []byte("="))
				if bytes.Equal(sub_info_part[0], []byte("AF")) {
					for _, tmp_af = range bytes.Split(sub_info_part[1], []byte(",")) {
						tmp_p, _ = strconv.ParseFloat(string(tmp_af), 32)
						af = append(af, float32(tmp_p))
					}
				}
				if bytes.Equal(sub_info_part[0], []byte("VT")) {
					vt = make([]byte, len(sub_info_part[1]))
					copy(vt, sub_info_part[1])
				}
			}

			var_prof_elem.AleFreq = append(var_prof_elem.AleFreq, 0)
			if len(af) == len(alt_arr) {
				alt_prob = float32(0)
				for i = 0; i < len(alt_arr); i++ {
					var_prof_elem.Variant = append(var_prof_elem.Variant, alt_arr[i])
					var_prof_elem.AleFreq = append(var_prof_elem.AleFreq, af[i])
					alt_prob += af[i]
				}
				var_prof_elem.AleFreq[0] = 1 - alt_prob
			} else {
				alt_prob = 1 / float32(1+len(alt_arr))
				for i = 0; i < len(alt_arr); i++ {
					var_prof_elem.Variant = append(var_prof_elem.Variant, alt_arr[i])
					var_prof_elem.AleFreq = append(var_prof_elem.AleFreq, alt_prob)
				}
				var_prof_elem.AleFreq[0] = alt_prob
			}
			var_prof_elem.VarType = vt
			pos = -1
			for i, tmp_name = range chr_name {
				if string(sub_line[0]) == string(tmp_name) {
					pos = chr_pos[i]
					break
				}
			}
			if pos == -1 {
				fmt.Println("Chromosome", string(sub_line[0]), "is missing from variant profile")
			}
			var_pos, _ = strconv.Atoi(string(sub_line[1]))
			var_prof[pos+var_pos-1] = var_prof_elem
		}
	}
	return var_prof
}
