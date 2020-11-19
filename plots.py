import matplotlib.pyplot as plt
import pandas as pd

file = open("C:\\Users\\Likhith\\source\\repos\\PDC_proj\\times.txt", "r")
a = file.read()

file2 = open("D:\\c programs\\PDC\\times.txt", "r")
b = file2.read()

file3 = open("D:\\c programs\\PDC2\\times.txt", "r")
c = file3.read()
# # analysing CUDA
# l = list(a.split("\n"))
# l1 = []
# file.close()
# for i in l:
#     l1.append(list(i.split()))

# df = pd.DataFrame(l1, columns = ['n', 'time', 'num'])
# df1 = df[6:12]
# num_matrices_3 = df1['num']
# time_3 = df1['time'].map(float)

# df2 = df[:6]
# num_matrices_4 = df2['num']
# time_4 = df2['time'].map(float)

# plt.plot(num_matrices_3, time_3, label = "dim 3")
# plt.plot(num_matrices_4, time_4, label = 'dim 4')
# plt.xlabel("No of matrices")
# plt.ylabel("Time taken to execute")
# plt.title("No of matrices VS time taken (CUDA)")
# plt.legend()
# plt.show()

# # print(df1)
# # print(df2)


# Analysing CPP code
# l = list(b.split("\n"))
# l1 = []
# file2.close()
# for i in l:
#     l1.append(list(i.split()))

# df = pd.DataFrame(l1, columns = ['n', 'time', 'num'])
# df1 = df[:6]
# num_matrices_3 = df1['num']
# time_3 = df1['time'].map(float)

# df2 = df[6:12]
# num_matrices_4 = df2['num']
# time_4 = df2['time'].map(float)

# plt.plot(num_matrices_3, time_3, label = "dim 3")
# plt.plot(num_matrices_4, time_4, label = 'dim 4')
# plt.xlabel("No of matrices")
# plt.ylabel("Time taken to execute")
# plt.title("No of matrices VS time taken (CPP)")
# plt.legend()
# plt.show()

#analyzing MPI code

# l_m = list(c.split("\n"))
# l1_m = []
# file3.close()
# for i in l_m:
#     l1_m.append(list(i.split()))

# df_m = pd.DataFrame(l1_m, columns = ['n', 'time', 'num'])
# df1_m = df_m[:6]
# num_matrices_3_m = df1_m['num']
# time_3_m = df1_m['time'].map(float)

# df2_m = df_m[6:12]
# num_matrices_4_m = df2_m['num']
# time_4_m = df2_m['time'].map(float)

# plt.plot(num_matrices_3_m, time_3_m, label = "dim 3")
# plt.plot(num_matrices_4_m, time_4_m, label = 'dim 4')
# plt.xlabel("No of matrices")
# plt.ylabel("Time taken to execute")
# plt.title("No of matrices VS time taken (MPI)")
# plt.legend()
# plt.show()



#comparing Both the codes
l = list(a.split("\n"))
l1 = []
file.close()
for i in l:
    l1.append(list(i.split()))
df = pd.DataFrame(l1, columns = ['n', 'time', 'num'])
df2 = df[:6]
num_matrices_4 = df2['num']
time_4 = df2['time'].map(float)


l2 = list(b.split("\n"))
l12 = []
file2.close()
for i in l2:
    l12.append(list(i.split()))
df12 = pd.DataFrame(l12, columns = ['n', 'time', 'num'])
df22 = df12[6:12]
num_matrices_42 = df22['num']
time_42 = df22['time'].map(float)

l_m = list(c.split("\n"))
l1_m = []
file3.close()
for i in l_m:
    l1_m.append(list(i.split()))

df_m = pd.DataFrame(l1_m, columns = ['n', 'time', 'num'])
df2_m = df_m[6:12]
num_matrices_4_m = df2_m['num']
time_4_m = df2_m['time'].map(float)

plt.plot(num_matrices_4, time_4, label = "CUDA")
plt.plot(num_matrices_42, time_42, label = 'CPP')
plt.plot(num_matrices_4_m, time_4_m, label = 'MPI')
plt.xlabel("No of matrices")
plt.ylabel("Time taken to execute")
plt.title("CUDA VS SERIAL-CPP vs MPI")
plt.legend()
plt.show()