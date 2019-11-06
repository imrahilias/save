import numpy as np
import re

A = np.array([[0, 0, 0, 1, 0, 5, 0, 6, 8],
              [0, 0, 0, 0, 0, 0, 7, 0, 1],
              [9, 0, 1, 0, 0, 0, 0, 3, 0],
              [0, 0, 7, 0, 2, 6, 0, 0, 0],
              [5, 0, 0, 0, 0, 0, 0, 0, 3],
              [0, 0, 0, 8, 7, 0, 4, 0, 0],
              [0, 3, 0, 0, 0, 0, 8, 0, 5],
              [1, 0, 5, 0, 0, 0, 0, 0, 0],
              [7, 9, 0, 4, 0, 1, 0, 0, 0]])

A=np.zeros((9,9))

sudoku=open("./sudokus/oe20060916_leicht.sudo","r")
Numbers=sudoku.readlines()
Numbers=re.sub("[^0-9]","",str(Numbers))
i=-1
for m in range(0,9):
    for n in range(0,9):
        i=i+1
        A[m,n]=Numbers[i]
print(A)







def get_list_of_possible_Numbers(Neighbours):
    list = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    list2 = [1, 2, 3, 4, 5, 6, 7, 8, 9]

    for i in list:
        if i in Neighbours:
            list2.remove(i)
    return list2


def get_Row(A,i,m):
    Row = A[:, m]
    return Row

def get_Column(A,i,m):
    Column = A[i, :]
    return Column

def get_Box(A, i, m):

    Newlist = []

    x = m
    y = i
    for k in range(0, 9):
        Newlist.append(A[y, x])
        if (x + 1) % 3 == 0:
            x = x - 2
            y = y + 1
            if (y) % 3 == 0:
                y = y - 3



        else:
            x = x + 1
    return Newlist


def find_one_Number(A):
    for m in range(0, 9):
        for n in range(0, 9):
            if A[m,n]==0:
                forbidden_Numbers = get_Box(A, m, n)
                forbidden_Numbers.extend(get_Row(A,m,n))
                forbidden_Numbers.extend(get_Column(A, m, n))
                print(forbidden_Numbers)
                possible_Numbers = get_list_of_possible_Numbers(forbidden_Numbers)
                if len(possible_Numbers) == 1:

                    print(A[m,n],forbidden_Numbers)
                    A[m, n] = possible_Numbers[0]
                    print("Hurray its a", (possible_Numbers[0]), "at", m, "and", n)
                else:
                    print("No Number found at [%s,%s]", (m,n))
                    print("Possible Numbers are: %s", forbidden_Numbers)

    return A

while 0 in A:
    A = find_one_Number(A)


print(A)
