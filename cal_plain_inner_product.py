#! /usr/bin/python

file1 = open("alice/p1.csv",'r')
file2 = open("bob/p2.csv",'r')

file1.readline()
file2.readline()

v1 = [float(ele.rstrip()) for ele in file1.readlines()]
v2 = [float(ele.rstrip()) for ele in file2.readlines()]

print  sum([v1[i]*v2[i] for i in range(len(v1))])

