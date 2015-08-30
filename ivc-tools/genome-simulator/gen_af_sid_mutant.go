/*----------------------------------------------------------------------------------
 Generate mutant genome based on SNPs and INDELs from VCF file
 By Nam S. Vo, March 2015
 Generate mutations (subs, indels) based on allele frequency
 Usage: go run generate_af_sid_mutant.go ref_file snp_prof_file result_dir
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
	"fmt"
	"os"
	"strconv"
	"bytes"
	"sort"
	"time"
	"path"
	"math/rand"
)

type SNPProfile struct {
	Profile    [][]byte
	AlleleFreq []float32
	SVType		string
}

func main(){
	if len(os.Args) < 4 {
		fmt.Println("Usage: go run genome-mutant.go <input.fasta> <input.vcf> <result_dir>")
		return
	}
	fmt.Println("Reading reference...")
	header, seq := ReadFASTA(os.Args[1])
	fmt.Println("Reading SNP profile...")
	SNP_Prof := ReadVCF(os.Args[2])
	fmt.Println("Generating mutant genome and its variant profile...")
	mutant, mut_var, all_mut_var := GenerateMutant(seq, SNP_Prof)
	fmt.Println("Saving mutant gennome and its variant profile...")
	SaveMutant(os.Args[3], header, mutant, mut_var, all_mut_var)
	fmt.Println("Done!")
}

//--------------------------------------------------------------------------------------------------
// Generate mutant genome
//--------------------------------------------------------------------------------------------------
func GenerateMutant(seq []byte, SNP_Prof map[int]SNPProfile) ([]byte, []map[int][][]byte, map[int][][]byte) {

	mutant := make([]byte, 0)
	mut_var := make([]map[int][][]byte, 6)
	for i := 0; i < 6; i++ {
		mut_var[i] = make(map[int][][]byte)
	}
	all_mut_var := make(map[int][][]byte)

	snp_pos := make([]int, 0, len(SNP_Prof))
	for key, _ := range SNP_Prof {
		snp_pos = append(snp_pos, key)

	}
	sort.Ints(snp_pos)
	fmt.Println("Total number of variants:", len(snp_pos))

	var pos, idx int

	rand.Seed(time.Now().UnixNano())
	mutant = append(mutant, seq[ : snp_pos[0]]...)

	var snp_prof [][]byte
	var allele_freq []float32
	var sv_type string
	sv_num, ol_var_num := 0, 0
	for idx = 0; idx < len(snp_pos) - 1; idx++ {
		pos = snp_pos[idx]
		snp_prof = SNP_Prof[pos].Profile
		allele_freq = SNP_Prof[pos].AlleleFreq
		sv_type = SNP_Prof[pos].SVType
		if sv_type != "SNP" && sv_type != "INDEL" { //ignore other types of variants
			mutant = append(mutant, seq[snp_pos[idx] : snp_pos[idx + 1]]...)
			sv_num++
		} else {
			if rand.Float32() < allele_freq[1] { //Take ALT bases
				if len(snp_prof[0]) == 1 { //SNPs or Insertions
					if len(snp_prof[1]) == 1 { //SNPs
						mut_var[0][pos] = snp_prof
					} else { //Insertions
						mut_var[1][pos] = snp_prof
					}
					all_mut_var[pos] = snp_prof
					mutant = append(mutant, snp_prof[1]...)
					mutant = append(mutant, seq[snp_pos[idx] + 1 : snp_pos[idx + 1]]...)
				} else if len(snp_prof[0]) > 1 { //Deletions
					mut_var[2][pos] = snp_prof
					all_mut_var[pos] = snp_prof
					mutant = append(mutant, snp_prof[1]...)
					//Ignore SNP profile that overlap with the selected deletion
					anchor_pos := snp_pos[idx] + len(snp_prof[0])
					for snp_pos[idx + 1] < anchor_pos {
						idx++
						ol_var_num++
					}
					mutant = append(mutant, seq[anchor_pos : snp_pos[idx + 1]]...)
				}
			} else { //Take REF bases
				if len(snp_prof[0]) == 1 { //SNPs or Insertions
					if len(snp_prof[1]) == 1 { //SNPs
						mut_var[3][pos] = append(mut_var[3][pos], snp_prof[0])
						mut_var[3][pos] = append(mut_var[3][pos], snp_prof[0])
					} else { //Insertions
						mut_var[4][pos] = append(mut_var[4][pos], snp_prof[0])
						mut_var[4][pos] = append(mut_var[4][pos], snp_prof[0])
					}
					all_mut_var[pos] = append(all_mut_var[pos], snp_prof[0])
					all_mut_var[pos] = append(all_mut_var[pos], snp_prof[0])
					mutant = append(mutant, snp_prof[0]...)
					mutant = append(mutant, seq[snp_pos[idx] + 1 : snp_pos[idx + 1]]...)
				} else if len(snp_prof[0]) > 1 { //Deletions
					mut_var[5][pos] = append(mut_var[5][pos], snp_prof[0])
					mut_var[5][pos] = append(mut_var[5][pos], snp_prof[0])
					all_mut_var[pos] = append(all_mut_var[pos], snp_prof[0])
					all_mut_var[pos] = append(all_mut_var[pos], snp_prof[0])
					mutant = append(mutant, snp_prof[0]...)
					//Ignore SNP profile that overlap with the selected deletion
					anchor_pos := snp_pos[idx] + len(snp_prof[0])
					for snp_pos[idx + 1] <  anchor_pos {
						idx++
						ol_var_num++
					}
					mutant = append(mutant, seq[anchor_pos : snp_pos[idx + 1]]...)
				}
			}
		}
	}
	mutant = append(mutant, seq[snp_pos[idx] : ]...) //ignore the last SNP at this moment
	fmt.Println("Total, sub_diff_ref, ins_diff_ref, del_diff_ref, sub_same_ref, ins_same_ref, del_same_ref, SV, OL_VAR:")
	fmt.Println(len(SNP_Prof), len(mut_var[0]), len(mut_var[1]), len(mut_var[2]), len(mut_var[3]), len(mut_var[4]), len(mut_var[5]), sv_num, ol_var_num)

	return mutant, mut_var, all_mut_var
}

//--------------------------------------------------------------------------------------------------
// Save mutant genome, mutation profiles to file
//--------------------------------------------------------------------------------------------------
func SaveMutant(dir_name string, header, mutant []byte, mut_var []map[int][][]byte, all_mut_var map[int][][]byte) {
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
	ref_file.Write(mutant)

	file_names := []string{"sub_diff_ref.txt", "ins_diff_ref.txt", "del_diff_ref.txt", "sub_same_ref.txt", "ins_same_ref.txt", "del_same_ref.txt"}
	for i, file_name := range file_names {
		WriteVarToFile(dir_name, file_name, mut_var[i])
	}
	WriteVarToFile(dir_name, "mut_var.txt", all_mut_var)
}

func WriteVarToFile(dir_name, file_name string, mut_var map[int][][]byte) {
	mut_var_file, err := os.Create(path.Join(dir_name, file_name))
	if  err != nil {
        return
    }
	defer mut_var_file.Close()

	var mut_var_key []int
	for key, _ := range mut_var {
		 mut_var_key = append(mut_var_key, key)
	}
	sort.Ints(mut_var_key)
	_, err = mut_var_file.WriteString("#POS\tREF\tALT\n");
	for _, pos := range mut_var_key {
		snp_pos := strconv.Itoa(pos)
		_, err = mut_var_file.WriteString(snp_pos + "\t" + string(mut_var[pos][0]) + "\t" + string(mut_var[pos][1]) + "\n")
		if err != nil {
            fmt.Println(err)
            break
        }
	}
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
	header := make([]byte, len(line))
	copy(header, line)
	var seq []byte
	for scanner.Scan() {
		line = scanner.Bytes()
		seq = append(seq, line...)
	}
	return header, seq
}

//--------------------------------------------------------------------------------------------------
// Read VCF files for human data
// 	Consider only biallelic mutations
// 	Take POS, REF, ALT, AF, and SVTYPE
//--------------------------------------------------------------------------------------------------
func ReadVCF(seq_file string) map[int]SNPProfile {
	SNP_Prof := make(map[int]SNPProfile)
	f, err := os.Open(seq_file)
	if err != nil {
		fmt.Printf("%v\n", err)
		os.Exit(1)
	}
	defer f.Close()

	var line, sub_info []byte
	var alt, info, sub_info_part [][]byte
	var pos int
	var af float32
	var tmp_af float64
	var tmp_svtype string

	rand.Seed(time.Now().UnixNano())
	data := bufio.NewReader(f)
	for {
		line, err = data.ReadBytes('\n')
		if err != nil {
			fmt.Println("Finishing reading file", seq_file)
			break
		}
		if bytes.Equal(line[0:1], []byte("#")) {
			//Do something here later
		} else {
			tmp := SNPProfile{}
			sub_line, _ := SplitN(line, []byte("\t"), 9)
			alt = bytes.Split(sub_line[4], []byte(","))
			if len(alt) > 1 {
				fmt.Println("Not a biallelic mutation", pos, len(alt))
				continue
			}
			info = bytes.Split(sub_line[7], []byte(";"))
			for _, sub_info = range info {
				sub_info_part = bytes.Split(sub_info, []byte("="))
				if bytes.Equal(sub_info_part[0], []byte("AF")) {
					tmp_af, _ = strconv.ParseFloat(string(sub_info_part[1]), 32)
					af = float32(tmp_af)
				}
				if bytes.Equal(sub_info_part[0], []byte("VT")) {
					tmp_svtype = string(sub_info_part[1])
				}
			}
			//Take REF seq
			tmp.Profile = append(tmp.Profile, sub_line[3])
			tmp.AlleleFreq = append(tmp.AlleleFreq, 1.0 - af)
			//Take ALT seq
			tmp.Profile = append(tmp.Profile, alt[0])
			tmp.AlleleFreq = append(tmp.AlleleFreq, af)
			//Take type of variants
			tmp.SVType = tmp_svtype
			//Take POS
			pos, _ = strconv.Atoi(string(sub_line[1]))
			//Record all info
			SNP_Prof[pos - 1] = tmp //convert 1-base index to 0-base index
		}
	}
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
