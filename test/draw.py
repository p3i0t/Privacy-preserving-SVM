#/usr/bin/python


import matplotlib.pyplot as plt
import pickle 


file = open('single_party_results.pickle','r')
encrypt_results = pickle.load(file)

x=[0]
enc_time=[0]
M_time=[0]
cal_time=[0]

for ele in encrypt_results:
	x.append(ele[0])
	time = [float(ele) for ele in ele[1].split()]
	enc_time.append(time[0])	
	M_time.append(time[1])	
	cal_time.append(time[2])	


plt.xlabel('Dimension of Plaintext')
plt.ylabel('Time comsumed/s')
#plt.title('Performance Comparison of Encryption')
plt.title('Performance of Optimized Encryption')
plt.plot(x,enc_time,marker='o',linestyle='--',color='r',label="encryption")
plt.plot(x,M_time,marker='^',linestyle='-.',color='blue',label="M generation")
plt.plot(x,cal_time,marker='+',linestyle='-.',color='blue',label="dot production")
plt.axis([0,500,0,10])
# add a legend, to let the viewer know which curve is which
# To do this, we should adding label arguments in plot method
# then call legend() to draw the legend box.
plt.legend(loc=9) 
#plt.show()
plt.savefig('Performance_of_single_party.png')


