#/usr/bin/python
# coding=utf-8

import subprocess
import pickle

def stat(arg):
	time_elapsed = 0.0
	count = 1
	for i in range(count):
		process = subprocess.Popen('./test '+str(arg), \
		shell= True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

		out, err = process.communicate()
		#time_elapsed +=  float(out)
	return out #time_elapsed/count


if __name__=='__main__':
	results = [(100*i,stat(100*i)) for i in range(1,6)]
	file = open('single_party_results.pickle','w')
	pickle.dump(results,file)
