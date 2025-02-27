import numpy as np
import scipy
import scipy.linalg
import time

def load_eigen_npy(filename):
    with open(filename, "rb") as f:
        rows = np.fromfile(f, dtype=np.int32, count=1)[0]
        cols = np.fromfile(f, dtype=np.int32, count=1)[0]
        data = np.fromfile(f, dtype=np.float64).reshape(rows, cols)
    return data

A = load_eigen_npy("/home/paul_johannssen/Desktop/masterarbeit/goast_old/build/examples/matrices/mult0.npy")
b = load_eigen_npy("/home/paul_johannssen/Desktop/masterarbeit/goast_old/build/examples/rhs/rhs0.npy")
sol = load_eigen_npy("/home/paul_johannssen/Desktop/masterarbeit/goast_old/build/examples/sol/sol0.npy")



start_time = time.time()
x = scipy.linalg.solve(A, b, assume_a="sym")
end_time = time.time()
# Compute the time taken
time_taken = end_time - start_time
print(f"time taken = {time_taken}")
print("Solution error python: ")
print(scipy.linalg.norm(A @ x - b))
print("Solution error c++: ")
print(scipy.linalg.norm(A @ sol - b))