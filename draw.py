#/usr/bin/python


import matplotlib.pyplot as plt
import pickle 


file = open('results','r')
encrypt_results = pickle.load(file)

print encrypt_results
x=range(0,51,10)
enc_time=[]
M_time=[]
inter_dec_time=[]
alice_dec_time=[]
bob_dec_time=[]
total_dec_time=[]

inner_time =[]
trans_time =[]
dec_time =[]

encrypt_results[0] = (0,0,0,0,0,0)
for ele in encrypt_results:
	#x.append(ele[0])
	time = list(ele)
	enc_time.append(time[0])
	
	M_time.append(time[1])	
	
	inter_dec_time.append(time[2])	
	alice_dec_time.append(time[3])	
	bob_dec_time.append(time[4])	
	total_dec_time.append(time[2]+time[3]+time[4])
	
	#inner_time.append(time[2])
	#trans_time.append(time[2])
	dec_time.append(time[5])
plt.xlabel('Dimension of Plaintext')
plt.ylabel('Time comsumed/s')
#plt.title('Performance Comparison of Encryption')
plt.title('Performance of System')
plt.plot(x,enc_time,marker='o',linestyle='--',color='r',label="encryption")

plt.plot(x,M_time,marker='^',linestyle='-.',color='blue',label="M generation")

plt.plot(x,inter_dec_time,marker='+',linestyle='-.',color='blue',label="inner production(inter)")
plt.plot(x,alice_dec_time,marker='1',linestyle='-.',color='black',label="inner production(alice)")
plt.plot(x,bob_dec_time,marker='*',linestyle='-.',color='green',label="inner production(bob)")
plt.plot(x,total_dec_time,marker='D',linestyle='-.',color='purple',label="inner production(total)")

#plt.plot(x,inner_time,marker='+',linestyle='-.',color='blue',label="inner product")
#plt.plot(x,trans_time,marker='D',linestyle='-.',color='black',label="transform")
plt.plot(x,dec_time,marker='*',linestyle='-.',color='green',label="decryption")

plt.axis([0,50,0,60])
# add a legend, to let the viewer know which curve is which
# To do this, we should adding label arguments in plot method
# then call legend() to draw the legend box.
plt.legend(loc=9) 
#plt.show()
plt.savefig('Performance_of_System.png')


