#/usr/bin/python
# coding=utf-8

import subprocess
import pickle
import os

data = [(0,0,0,0,0)]

for i in range(100,501,100):
	enc_time = 0.0
	M_time = 0.0
	inter_dec_time = 0.0
	alice_dec_time = 0.0
	bob_dec_time = 0.0
	inner_time = 0.0
	trans_time = 0.0
	dec_time =0.0

	os.chdir('alice')
	process = subprocess.Popen('./alice '+str(10)+" "+str(i), \
	shell= True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out,err = process.communicate()
	time = [float(ele)for ele in out.split()]
	enc_time += time[0]
	M_time += time[1] + time[2]

	os.chdir('../bob')
	process = subprocess.Popen('./bob '+str(10)+" "+str(i), \
	shell= True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out,err = process.communicate()
	time = [float(ele)for ele in out.split()]
	enc_time += time[0]
	M_time += time[1] + time[2]

	os.chdir('../cloud')
	process = subprocess.Popen('./cloud '+str(3)+" "+str(10), \
	shell= True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out,err = process.communicate()
	time = [float(ele)for ele in out.split()]
	M_time += time[0]
	'''
	inter_dec_time += time[1]
	alice_dec_time += time[2]
	bob_dec_time += time[3]
	print i, " cnt: ", time[4], time[5], time[6]
	'''
	inner_time += time[1]
	trans_time += time[2]
	dec_time += time[3]

	os.chdir('..')
	#data.append((enc_time,M_time,inter_dec_time,alice_dec_time,bob_dec_time))
	data.append((enc_time,inner_time,trans_time,dec_time))
print data
file = open("results",'w')
pickle.dump(data,file)
