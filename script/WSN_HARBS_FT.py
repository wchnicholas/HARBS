#!/usr/bin/python
import numpy as np
from scipy.linalg import hadamard

def main():
  M=3; n=2**M
  H=hadamard(n)
  print float(1)/float(8)*np.dot(H,np.array([[73.9],[77.5],[71.2],[74.5],[73.0],[75.5],[70.6],[73.4]]))

if __name__ == "__main__":
  main()
