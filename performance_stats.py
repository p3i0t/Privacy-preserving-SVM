#/usr/bin/python
# coding=utf-8

import subprocess
import pickle
import os

data = [(0,0,0,0,0)]

for i in range(10,11,10):
	enc_time = 0.0
	M_time = 0.0
	inter_dec_time = 0.0
	alice_dec_time = 0.0
	bob_dec_time = 0.0
	inner_time = 0.0
	trans_time = 0.0
	dec_time =0.0

	os.chdir('alice')
	process = subprocess.Popen('./alice '+str(i)+" "+str(200), \
	shell= True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out,err = process.communicate()
	time = [float(ele)for ele in out.split()]
	enc_time += time[0]
	M_time += time[1] + time[2]

	os.chdir('../bob')
	process = subprocess.Popen('./bob '+str(i)+" "+str(200), \
	shell= True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out,err = process.communicate()
	time = [float(ele)for ele in out.split()]
	enc_time += time[0]
	M_time += time[1] + time[2]

	os.chdir('../cloud')
	process = subprocess.Popen('./cloud '+str(3)+" "+str(i), \
	shell= True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out,err = process.communicate()
	time = [float(ele)for ele in out.split()]
	M_time += time[0]
	
	inter_dec_time += time[1]
	dec_time += time[2]
	alice_dec_time += time[3]
	dec_time += time[4]
	bob_dec_time += time[5]
	dec_time += time[6]
	#print i, " cnt: ", time[4], time[5], time[6]
	'''
	inner_time += time[1]
	trans_time += time[2]
	dec_time += time[3]
	'''
	os.chdir('..')
	#data.append((enc_time,M_time,inter_dec_time,alice_dec_time,bob_dec_time))
	data.append((enc_time,M_time,inter_dec_time,alice_dec_time,bob_dec_time,dec_time))
print data
file = open("results",'w')
pickle.dump(data,file)
