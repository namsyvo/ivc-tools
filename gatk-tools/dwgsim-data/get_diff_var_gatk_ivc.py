"""
Get alignment info for fp (snp, indel)
Input: result folder, length of read/ref/mut
Ouput: alignment info (read-ref)
Usage: python get_diff_var_gatk_ivc.py cov_num var_type ivc_call_var_dir
"""
import sys
import os

if __name__ == "__main__":

    if len(sys.argv) != 4:
        print "Usage: python get_diff_var_gatk_ivc.py cov_num var_type ivc_call_var_dir"
        exit(0)

    cov_num = int(sys.argv[1])
    var_type = int(sys.argv[2])
    ivc_result_path = sys.argv[3]

    ref_len = 249250621
    read_lens = 100
    read_nums = str(cov_num*ref_len/(2*read_lens))

    '''
    Find variants called by GATK but not by IVC
    '''
    fp_fn = ["dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_snp_unknown.1.14.txt", \
             "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_snp_known.1.14.txt", \
             "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_indel_unknown.1.14.txt", \
             "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_indel_known.1.14.txt"]

    ivc_result_dn = os.path.join("/data/nsvo/test-data/GRCh37_chr1/results/sim-reads/af_sid_mutant_dwgsim/ivc_0.70", ivc_result_path, "fpfntp_info")
    var_pos = {}
    fp_inf = open(os.path.join(ivc_result_dn, fp_fn[var_type]))
    line = fp_inf.readline()
    for line in fp_inf:
        tmp = line.strip().split()
        var_pos[tmp[0]] = True

    fp_fn = ["dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_snp_unknown.20.0.txt", \
             "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_snp_known.20.0.txt", \
             "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_indel_unknown.20.0.txt", \
             "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_indel_known.20.0.txt"]

    gatk_result_dn = "/data/nsvo/test-data/GRCh37_chr1/results/sim-reads/af_sid_mutant_dwgsim/gatk/fpfntp_info"
    diff_var_outf = open(os.path.join(gatk_result_dn, fp_fn[var_type] + "-diff-var-gatk-not-ivc"), "w")
    fp_inf = open(os.path.join(gatk_result_dn, fp_fn[var_type]))
    line = fp_inf.readline()
    for line in fp_inf:
        tmp = line.strip().split()
        if tmp[0] not in var_pos:
            diff_var_outf.write(line)
    diff_var_outf.close()

    '''
    Find variants called by IVC but not by GATK
    '''
    fp_fn = ["dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_snp_unknown.20.0.txt", \
             "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_snp_known.20.0.txt", \
             "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_indel_unknown.20.0.txt", \
             "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_indel_known.20.0.txt"]

    gatk_result_dn = "/data/nsvo/test-data/GRCh37_chr1/results/sim-reads/af_sid_mutant_dwgsim/gatk/fpfntp_info"
    var_pos = {}
    fp_inf = open(os.path.join(gatk_result_dn, fp_fn[var_type]))
    line = fp_inf.readline()
    for line in fp_inf:
        tmp = line.strip().split()
        var_pos[tmp[0]] = True

    fp_fn = ["dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_snp_unknown.1.14.txt", \
             "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_snp_known.1.14.txt", \
             "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_indel_unknown.1.14.txt", \
             "dwgsim_reads_100.0.00015-0.0015." + read_nums + ".tp_indel_known.1.14.txt"]

    ivc_result_dn = os.path.join("/data/nsvo/test-data/GRCh37_chr1/results/sim-reads/af_sid_mutant_dwgsim/ivc_0.70", ivc_result_path, "fpfntp_info")
    diff_var_outf = open(os.path.join(ivc_result_dn, fp_fn[var_type] + "-diff-var-ivc-not-gatk"), "w")
    diff_var_unique_outf = open(os.path.join(ivc_result_dn, fp_fn[var_type] + "-diff-var-unique_ivc-not-gatk"), "w")
    fp_inf = open(os.path.join(ivc_result_dn, fp_fn[var_type]))
    line = fp_inf.readline()
    ivc_var_pos = {}
    prev_pos = ""
    for line in fp_inf:
        tmp = line.strip().split()
        if tmp[0] not in var_pos:
            ivc_var_pos[tmp[0]] = True
            diff_var_outf.write(line)
            if tmp[0] != prev_pos:
                diff_var_unique_outf.write(line)
            prev_pos = tmp[0]
    diff_var_outf.close()
    diff_var_unique_outf.close()
