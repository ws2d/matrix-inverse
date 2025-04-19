# Efficient Algorithms for Computing the Inverse of Encrypted Matrices

## Compile and Run (Ubuntu 20.04.6 LTS)

### Dependencies

\- cmake 3.13 or higher

\- clang

\- [Microsoft SEAL](https://github.com/microsoft/seal)

### Compile

Download the source code from GitHub, unzip and enter the repository root directory, and then run the following:

  cd example

  mkdir build

  cd build

  cmake ../test/

  make

###  How to use

\- Open a terminal, enter the `build` directory and run 

   ./test_example

### Test Instructions

Once the file is executed, it will output the following information:

1. triangle inv(4,4)

2. dia: pack new encoding(d,d)

3. real: pack new encoding(d,d)

4. tae inv(d,d)

5. dia:block inversion(256,256)

6. dia:block inversion(512,512)

7--10. test

0. Exit



In this case, `case 1-6` refer to the experiment, while `case 7-10` are our tests. The introduction of our files are as follows:

It should be pointed out that, unless otherwise specified, the files under the build folder are encoded by entry by entry, and the files under test are encoded by row by row.



build:

All results is achieved in  encryped domain.

dia_block and dia_result:  The result of  diagonal dominant matrices of  different dimension and (256,256),(512,512) is row by row.

dia_input: The input of diagonal dominant matrices of different dimension

dia_movie:The input and result of ml-100k dataset in the inverse of E step.

dia_triangle:  The input and result of three methods in  Table 4.(triangle method is row by row)

real_input: The input of real symmetric matrices of different dimension.

real_movie: The input and result of ml-100k dataset in the inverse of M step.

real_result: The result of real symmetric matrices of  different dimension.

real_tae_A:The  input of tae method(ALYY), corresponding to "real_input" with $\alpha$.

real_tae_result: The result of ALYY.

real_tae_Y0: The Y0 of ALYY.



src: 

Source code files.



test: 

dia_A: The matrix A of diagonal dominant matrices, corresponding to "dia_input" with  $\alpha$.

dia_alpha: The $\alpha$ of diagonal dominant matrices.

dia_input: The input of diagonal dominant matrices of different dimension.

dia movie: The input of ml-100k dataset in the inverse of E step.

dia_triangle:  The input and result of three methods in  Table 4.

real_A: The matrix A of real symmetric matrices, corresponding to "real_input" with  $\alpha$.

real_alpha: The $\alpha$ of real symmetric matrices.

real_input: The input of real symmetric matrices of different dimension.

real_movie: The input of ml-100k dataset in the inverse of M step.

real_tae_Y0: The Y0 of ALYY.

test_example.cpp: The main function.



### Contact Me

Because different experiments may use the same function, there will be slight adjustments to the input and output of the function. If you have any questions, please contact me.

Email:  622230830030@mails.cqjtu.edu.cn



















