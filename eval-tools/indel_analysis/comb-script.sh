python get_tpfn_htd_indels.py ../../ivc-tools/0.8.0-test-scripts/config_chrall.json 0.70 tp_snp /data/nsvo/test_data/GRCh37/results/sim_reads/af_sid_mutant_dwgsim/ivc_0.70/IVC-0.8.1/fpfntp_info/dwgsim_reads_100.0.00015-0.0015.464351611.tp_snp_unknown.20.0.10.0.txt ivc_tp_snp_htd.txt &
python get_tpfn_htd_indels.py ../../ivc-tools/0.8.0-test-scripts/config_chrall.json 0.70 tp_indel /data/nsvo/test_data/GRCh37/results/sim_reads/af_sid_mutant_dwgsim/ivc_0.70/IVC-0.8.1/fpfntp_info/dwgsim_reads_100.0.00015-0.0015.464351611.tp_indel_unknown.20.0.10.0.txt ivc_tp_indel_htd.txt &
python get_tpfn_htd_indels.py ../../ivc-tools/0.8.0-test-scripts/config_chrall.json 0.70 fn_snp /data/nsvo/test_data/GRCh37/results/sim_reads/af_sid_mutant_dwgsim/ivc_0.70/IVC-0.8.1/fpfntp_info/dwgsim_reads_100.0.00015-0.0015.464351611.fn_snp_unknown.20.0.10.0.txt ivc_fn_snp_htd.txt &
python get_tpfn_htd_indels.py ../../ivc-tools/0.8.0-test-scripts/config_chrall.json 0.70 fn_indel /data/nsvo/test_data/GRCh37/results/sim_reads/af_sid_mutant_dwgsim/ivc_0.70/IVC-0.8.1/fpfntp_info/dwgsim_reads_100.0.00015-0.0015.464351611.fn_indel_unknown.20.0.10.0.txt ivc_fn_indel_htd.txt &


python get_tpfn_htd_indels.py ../../ivc-tools/0.8.0-test-scripts/config_chrall.json 0.70 tp_snp /data/nsvo/test_data/GRCh37/results/sim_reads/af_sid_mutant_dwgsim/gatk_hc_realign/fpfntp_info/dwgsim_reads_100.0.00015-0.0015.464351611.tp_snp_unknown.20.0.txt hc_tp_snp_htd.txt &
python get_tpfn_htd_indels.py ../../ivc-tools/0.8.0-test-scripts/config_chrall.json 0.70 tp_indel /data/nsvo/test_data/GRCh37/results/sim_reads/af_sid_mutant_dwgsim/gatk_hc_realign/fpfntp_info/dwgsim_reads_100.0.00015-0.0015.464351611.tp_indel_unknown.20.0.txt hc_tp_indel_htd.txt &
python get_tpfn_htd_indels.py ../../ivc-tools/0.8.0-test-scripts/config_chrall.json 0.70 fn_snp /data/nsvo/test_data/GRCh37/results/sim_reads/af_sid_mutant_dwgsim/gatk_hc_realign/fpfntp_info/dwgsim_reads_100.0.00015-0.0015.464351611.fn_snp_unknown.20.0.txt hc_fn_snp_htd.txt &
python get_tpfn_htd_indels.py ../../ivc-tools/0.8.0-test-scripts/config_chrall.json 0.70 fn_indel /data/nsvo/test_data/GRCh37/results/sim_reads/af_sid_mutant_dwgsim/gatk_hc_realign/fpfntp_info/dwgsim_reads_100.0.00015-0.0015.464351611.fn_indel_unknown.20.0.txt hc_fn_indel_htd.txt &


cat ivc_tp_snp_htd.txt.K ivc_tp_indel_htd.txt.K ivc_fn_snp_htd.txt.K ivc_fn_indel_htd.txt.K > ivc_tpfn_htd.txt.K
cat ivc_tp_snp_htd.txt.U ivc_tp_indel_htd.txt.U ivc_fn_snp_htd.txt.U ivc_fn_indel_htd.txt.U > ivc_tpfn_htd.txt.U


cat hc_tp_snp_htd.txt.K hc_tp_indel_htd.txt.K hc_fn_snp_htd.txt.K hc_fn_indel_htd.txt.K > hc_tpfn_htd.txt.K
cat hc_tp_snp_htd.txt.U hc_tp_indel_htd.txt.U hc_fn_snp_htd.txt.U hc_fn_indel_htd.txt.U > hc_tpfn_htd.txt.U
