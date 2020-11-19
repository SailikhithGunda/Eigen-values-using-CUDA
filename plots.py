import matplotlib.pyplot as plt
import pandas as pd

file = open("C:\\Users\\Likhith\\source\\repos\\PDC_proj\\times.txt", "r")
a = file.read()

# analysing CUDA
l = list(a.split("\n"))
l1 = []
file.close()
for i in l:
    l1.append(list(i.split()))

df = pd.DataFrame(l1, columns = ['n', 'time', 'num'])
df1 = df[6:12]
num_matrices_3 = df1['num']
time_3 = df1['time'].map(float)

df2 = df[:6]
num_matrices_4 = df2['num']
time_4 = df2['time'].map(float)

plt.plot(num_matrices_3, time_3, label = "dim 3")
plt.plot(num_matrices_4, time_4, label = 'dim 4')
plt.xlabel("No of matrices")
plt.ylabel("Time taken to execute")
plt.title("No of matrices VS time taken (CUDA)")
plt.legend()
plt.show()

