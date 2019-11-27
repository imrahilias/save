#! /usr/bin/python
#-*- coding: utf-8 -*-

# init
import os
cpu = [0]
print("current, diff")

for i in range(1,9):
    cmd = os.popen("cat /proc/stat | grep 'cpu ' | awk '{print $2}\'").read()
    cputime = int(cmd);
    cpu.append(cputime)
    diff = cpu[i] - cpu[i-1]
    print(cputime, diff)
