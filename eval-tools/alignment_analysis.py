result_path = "/home/genomes/SNP_calling/results/soy_bean/chr1" + "/2-28-2014"
f=open(result_path + "/test100k-prog_stat-100.0.01.3x.txt")
data=f.readlines()
print len(data)
left_m_arr, left_n_arr, right_m_arr, right_n_arr = [], [], [], []

left_m_count, left_n_count, right_m_count, right_n_count = 0, 0, 0, 0
left_diff_count, right_diff_count, left_same_count, right_same_count = 0, 0, 0, 0

for line in data:
	nums = line.split()
	#print nums
	left_m_arr.append(int(nums[0]))
	left_n_arr.append(int(nums[1]))
	right_m_arr.append(int(nums[2]))
	right_n_arr.append(int(nums[3]))
	if int(nums[0]) == 0:
		left_m_count += 1
	if int(nums[1]) == 0:
		left_n_count += 1
	if int(nums[2]) == 0:
		right_m_count += 1
	if int(nums[3]) == 0:
		right_n_count += 1

	if int(nums[0]) != int(nums[1]):
		left_diff_count += 1
	if int(nums[2]) != int(nums[3]):
		right_diff_count += 1
	if int(nums[0]) == int(nums[1]):
		left_same_count += 1
	if int(nums[2]) == int(nums[3]):
		right_same_count += 1

print "left_m_count", left_m_count
print "left_n_count", left_n_count
print "right_m_count", right_m_count
print "right_n_count", right_n_count

print left_diff_count
print right_diff_count
print left_same_count - left_m_count
print right_same_count - right_m_count
