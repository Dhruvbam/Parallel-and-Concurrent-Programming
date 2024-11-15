# Parallel and Concurrent Programming
![image](https://github.com/Dhruvbam/Parallel-and-Concurrent-Programming/blob/main/Images/ss.png)

## About
This repository contains five assignment projects completed as part of the Parallel and Concurrent Programming (CS 4379) course at Texas Tech University. These projects focus on applying parallel computing techniques using **MPI (Message Passing Interface)** to solve computational problems efficiently. The assignments explore topics like parallel computation, matrix operations, synchronization, and force calculations, demonstrating various parallel programming paradigms.

## Assignments Overview

1. **Assignment 1: Parallel Row Sum of a Matrix**
   - **Topic:** MPI and Parallel Matrix Operations
   - **Description:** This project involves parallelizing the computation of row sums for a 100x100 matrix using two processes. The work is divided between the processes to reduce computation time by overlapping communication with computation.
   - [View Code](https://github.com/Dhruvbam/Parallel-and-Concurrent-Programming/blob/main/Assignment%201%20-%20Parallel%20Row%20Sum%20of%20a%20Matrix/Assignment_1.c)

2. **Assignment 2: Parallelizing Nested Loops**
   - **Topic:** MPI and Parallel Loop Distribution
   - **Description:** This project focuses on splitting the iteration space of nested loops across multiple processes using MPI to parallelize the workload. Processes perform independent computations and communicate to find the global least value.
   - [View Code](https://github.com/Dhruvbam/Parallel-and-Concurrent-Programming/blob/main/Assignment%202%20-%20Parallelizing%20Nested%20Loops/Assignment_2.c)

3. **Assignment 3: Custom Barrier Synchronization**
   - **Topic:** MPI and Synchronization
   - **Description:** Implements a custom barrier synchronization mechanism using MPI to ensure that all processes reach a synchronization point before continuing. The project emphasizes minimizing synchronization steps through efficient message passing.
   - [View Code](https://github.com/Dhruvbam/Parallel-and-Concurrent-Programming/blob/main/Assignment%203%20-%20Custom%20Barrier%20Synchronization/Assignment_3.c)

4. **Assignment 4: Parallel Matrix Multiplication**
   - **Topic:** MPI and Grid-Based Matrix Operations
   - **Description:** This project parallelizes a matrix multiplication algorithm using a grid-based decomposition where submatrices are assigned to different processes. Communication between processes ensures that row and column data are shared across iterations.
   - [View Code](https://github.com/Dhruvbam/Parallel-and-Concurrent-Programming/blob/main/Assignment%204%20-%20Parallel%20Matrix%20Multiplication/Assignment_4.c)

5. **Assignment 5: Parallel Force Calculation**
   - **Topic:** MPI and Load Balancing
   - **Description:** This project implements a parallel algorithm to calculate net forces on particles in space. The workload is balanced across processes by dividing the triangular-shaped work into equal sections, using log₂(p) communication complexity.
   - [View Code](https://github.com/Dhruvbam/Parallel-and-Concurrent-Programming/blob/main/Assignment%205%20-%20Parallel%20Force%20Calculation/Assignment_5.c)


## Built With
This repository primarily utilizes:
- <a href="https://www.open-mpi.org/" target="_blank" rel="noreferrer"><img src="https://img.shields.io/badge/MPI-00599C?style=for-the-badge&logo=mpi&logoColor=white" width="36" height="36" alt="MPI" /></a> **MPI (Message Passing Interface)**: Used for parallelization and inter-process communication.
- <a href="https://en.wikipedia.org/wiki/C_(programming_language)" target="_blank" rel="noreferrer"><img src="https://img.shields.io/badge/C-00599C?style=for-the-badge&logo=c&logoColor=white" width="36" height="36" alt="C Programming" /></a> **C**: The core programming language used to implement parallel algorithms.

## How to Use
1. Clone the repository:
    ```bash
    git clone https://github.com/your-repo/parallel-and-concurrent-programming.git
    ```
2. Navigate to the assignment folder of your choice.
3. Compile and run the C program using `mpicc` (MPI compiler):
    ```bash
    mpicc Assignment_1.c -o Assignment_1
    mpirun -np 2 ./Assignment_1
    ```

## Learning Outcomes
Through these assignments, I developed strong skills in parallel computing and system synchronization, focusing on:

- **Parallel Algorithms with MPI**: Implemented parallel algorithms using **MPI**, optimizing performance and resource utilization across distributed systems.
- **Efficient Synchronization**: Designed and implemented synchronization mechanisms to ensure data consistency across processes, mastering mutexes and semaphores.
- **Workload Balancing**: Gained expertise in balancing workloads across multiple processes, improving system efficiency in parallel computations.
- **Complex Communication Patterns**: Handled complex inter-process communication, enhancing my ability to design scalable and efficient distributed systems.

These experiences solidified my ability to tackle parallel computing challenges in real-world environments.

