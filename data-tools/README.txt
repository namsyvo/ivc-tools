#Verify ref and dbsnp files
python verify_ref.py
python verify_ref_dbsnp.py

#Generate reads
python generate_reads.py
python verify_reads.py
python convert_to_fq.py (will call convert_to_fq.go)
