package main

import "fmt"

func main() {
	//				0	1	 2	  3	   4    5    6    7    8    9    10
	//				   'C', 'A', 'T', 'A', 'C'
	seq := []byte{'A', 'C', 'A', 'T', 'A', 'C', 'G', 'C', 'A', 'T', 'T'}
	snp_pos := []int{1, 2, 3, 7, 9}
	mutant := make([]byte, 0)
	idx := 0
	mutant = append(mutant, seq[ : snp_pos[idx]]...)
	for idx=0; idx<1; idx++ {
		snp0 := []byte{'C', 'A', 'T', 'A', 'C'}
		//snp1 := []byte{'C'}
		mutant = append(mutant, snp0...)
		//Ignore SNP profile that overlap with the selected deletion
		//anchor_idx := idx
		anchor_pos := snp_pos[idx] + len(snp0)
		fmt.Println("anchor_pos", anchor_pos)
		i_d := 0
		for snp_pos[idx + 1] < anchor_pos {
			i_d++
			idx++
		}
		fmt.Println("next idx", idx, snp_pos[idx], string(seq[anchor_pos : snp_pos[idx + 1]]))
		fmt.Println("#ignored_indels", i_d)
		mutant = append(mutant, seq[anchor_pos : snp_pos[idx + 1]]...)
	}
	mutant = append(mutant, seq[snp_pos[idx] : ]...)
	fmt.Println(string(mutant))
}
