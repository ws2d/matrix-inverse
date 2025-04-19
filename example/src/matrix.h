#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <seal/seal.h>
#include <fstream>
#include <numeric>
#include <seal/util/numth.h>


using namespace seal::util;

template <typename T>
class Matrix
{
protected:
    std::size_t n, d;
    std::vector<std::vector<T>> M;

public:

    // matrix();

    // empty matrix
    Matrix() : n(0), d(0) {}

    //~matrix();

    //  rows * cols, all elements are initialized with the default constructor of T
    Matrix(std::size_t rows, std::size_t cols) : n(0), d(0)
    {
        resize(rows, cols);
    }

    void clear()
    {
        n = d = 0;
        M.clear();
    }

    void readFromFile(const std::string& filename)
    {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        for (std::size_t i = 0; i < n; i++)
        {
            for (std::size_t j = 0; j < d; j++)
            {
                T value;
                if (!(file >> value)) {
                    std::cerr << "Error reading value at position (" << i << ", " << j << ")" << std::endl;
                    file.close();
                    return;
                }
                set(j,i, value); // 将读取到的值存储到矩阵中
            }
        }

        file.close();
    }

    void resize(std::size_t rows, std::size_t cols)
    {
        std::size_t j;
        n = rows;
        d = cols;
        M.resize(d);
        for (j = 0; j < d; j++)
        {
            M[j].resize(n);
        }
    }

    // return the number of rows
    size_t get_rows() const
    {
        return n;
    }

    // return the number of columns
    size_t get_cols() const
    {
        return d;
    }


    // a reference to the element (i, j)
    T &operator()(const std::size_t i, const std::size_t j) { return M[j][i]; }

    Matrix<T> add(Matrix<T> & N){
        std::size_t m = N.get_rows();
        std::size_t p = N.get_cols();
        Matrix<T> C;
        C.resize(n, d);
        if ((m != n) || (d != p)){
            std::cerr << "matrix::add: the dimensions do not match. " << std::endl;
        }
        else{
            for (std::size_t i=0; i < n; i++){
                for (std::size_t j = 0; j<d; j++){
                    C.set(i, j, M[j][i] + N(i, j));
                }
            }
        }
        return C;
    }

    Matrix<T> multiply(Matrix<T> & N){
        std::size_t m = N.get_rows();
        std::size_t p = N.get_cols();
        Matrix<T> C;
        C.resize(n, p);
        if (m != d){
            std::cerr << "matrix::multiply: the dimensions do not match. " << std::endl;
        }
        else{
            for (std::size_t i=0; i < n; i++){
                for (std::size_t j = 0; j<p; j++){
                    T t = 0;
                    for(std::size_t k = 0; k < m; k++){
                        t += M[k][i] * N(k, j);
                    }
                    C.set(i, j, t);
                }
            }
        }
        return C;
    }

    inline void set(const std::size_t i, const std::size_t j, const T a) { M[j][i] = a; }
    
    inline T &get(const std::size_t i, const std::size_t j) { return M[j][i]; }

    inline T &gets(const std::size_t i, const std::size_t j) { return M[i][j]; }
    
    void Merage_matrix(size_t start_row,size_t end_row,size_t start_col,size_t end_col,Matrix<T> B){
        for(size_t i=start_row;i<end_row;i++){
            for(size_t j=start_col;j<end_col;j++){
                M[j][i]=B.get(i-start_row,j-start_col);
            }
        }
    }
    
    void print_pack(std::size_t rows, std::size_t cols)
    {
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                std::cout << std::setw(6) << std::right << std::right << M[i][j] << ",";
            }
        }
    }

    void print_1(std::size_t rows = 6, std::size_t cols = 6)
    {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < d; j++) {
                std::cout << std::setw(6) << std::right << std::right << M[i][j] << ",";
            }
        }
    }
    void print(std::size_t rows = 6, std::size_t cols = 6)
    {
        if ((n > 2*rows) && (d > 2*cols)){
            size_t r = rows / 2;
            size_t c = cols / 2;

            for (size_t i = 0; i < r; i++){
                std::cout << "    [";
                for (size_t j = 0; j < c; j++){
                    std::cout << std::setw(6) << std::right << std::right << M[j][i] << ",";
                }
                std::cout << std::setw(6) << std::right << std::right << "..."
                      << ",";
                for (size_t j = d - cols + c; j < d - 1; j++){
                    std::cout << std::setw(6) << std::right << std::right << M[j][i] << ",";
                }
                std::cout << std::setw(6) << std::right << std::right << M[d - 1][i] << "]" << std::endl;
            }

            std::cout << "    [";
            for (size_t j = 0; j < c; j++){
                std::cout << std::right << std::setw(6) << std::right << "..."
                      << ",";
            }
            std::cout << std::setw(6) << std::right << std::right << "..."
                  << ",";
            for (size_t j = d - cols + c; j < d - 1; j++){
                std::cout << std::setw(6) << std::right << std::right << "..."
                      << ",";
            }
            std::cout << std::setw(6) << std::right << std::right << "..."
                  << "]" << std::endl;

            for (size_t i = n - rows + r; i < n; i++){
                std::cout << "    [";
                for (size_t j = 0; j < c; j++){
                    std::cout << std::setw(6) << std::right << std::right << M[j][i] << ",";
                }
                std::cout << std::setw(6) << std::right << std::right << "..."
                      << ",";
                for (size_t j = d - cols + c; j < d - 1; j++){
                    std::cout << std::setw(6) << std::right << std::right << M[j][i] << ",";
                }
                std::cout << std::setw(6) << std::right << std::right << M[d - 1][i] << "]" << std::endl;
            }
        }
        else{
            for (size_t i = 0; i < n; i++){
                std::cout << "    [";
                for (size_t j = 0; j < d-1; j++){
                    std::cout << std::setw(6) << std::right << std::right << M[j][i] << ",";
                }
                std::cout << std::setw(6) << std::right << std::right << M[d - 1][i] << "]" << std::endl;
            }
        }
    }

    // flatten matrix to a vector by rows 
    std::vector<T> flatten_matrix_to_rows_vector() {
        std::vector<T> flat_vector(n*d);

        for (std::size_t i = 0; i < n; i++) {
            for (std::size_t j = 0; j < d; j++) {
                flat_vector[i*d+j] = M[j][i] ;
            }
        }
        return flat_vector;
    }

    //Jiang function generate matrix
    //generate u_sigma matrix
    void generate_u_sigma(std::size_t rows, std::size_t cols) {
        resize(rows * cols, rows * cols);
        for (std::size_t i = 0; i < cols; i++) {
            for (std::size_t j= 0; j < rows; j++) {
                M[i * rows + (j + i) % rows][i * rows + j] = 1;
            }
        }
    }

    //generate u_tau matrix
    void generate_u_tau(std::size_t rows, std::size_t cols) {
        resize(rows * cols, rows * cols);
        for (std::size_t i = 0; i < cols; i++) {
            for (std::size_t j= 0; j < rows; j++) {
                M[((rows + 1) * j + cols * i) % (rows * cols)][i * rows + j] = 1;
                // std::cout<<"rows:"<<i * rows + j<<"  cols:"<<((rows + 1) * j + cols * i) % (rows * cols)<<std::endl;
            }
        }
    }

    // get i-th diag vector
    std::vector<T> diag_vector(std::size_t i,size_t slot_conunt) {
        std::vector<T> diag_vector;
        if (n != d) {
            std::cout << "matrix is not a squre" << std::endl;
            return diag_vector;
        }
        for (std::size_t j = 0; j < n; j++) {
            diag_vector.push_back(M[(i + j) % d][j]);
        }
        diag_vector.resize(slot_conunt);
        return diag_vector;
    }

    void scaling_inplace(T c){
        for (std::size_t i=0; i < n; i++){
            for (std::size_t j = 0; j<d; j++){
                M[j][i] *= c;
            }
        }
    }

    bool is_zero(int64_t a = -5){
        double large_error=0;
        for (std::size_t i=0; i < n; i++){
            for (std::size_t j = 0; j<d; j++){
                if (abs(double((M[j][i]))) > (double)pow(10, a)){
                    std::cout << "Largest absolute error > "<<pow(10,a) << std::endl;
                    return false;
                }
                if(abs(double((M[j][i])))>large_error){
                    large_error=abs(double((M[j][i])));
                }
            }
        }
        std::cout<<"Largest absolute error: |"<<std::setprecision(10)<<large_error<<"|"<<std::endl;
        return true;
    }

    // 分割矩阵为四个子矩阵
    void split(size_t d,Matrix<T>& subMatrix1, Matrix<T>& subMatrix2, Matrix<T>& subMatrix3, Matrix<T>& subMatrix4) const {
        std::size_t subSize = d;

        for (std::size_t i = 0; i < subSize; ++i) {
            for (std::size_t j = 0; j < subSize; ++j) {
                subMatrix1.set(j, i, M[i][j]);
                subMatrix2.set(j, i, M[i][j + subSize]);
                subMatrix3.set(j, i, M[i + subSize][j]);
                subMatrix4.set(j, i, M[i + subSize][j + subSize]);
            }
        }
    }

    // 分割矩阵为十六个子矩阵
    void split(size_t d,Matrix<T>& subMatrix1, Matrix<T>& subMatrix2, Matrix<T>& subMatrix3, Matrix<T>& subMatrix4,
    Matrix<T>& subMatrix5, Matrix<T>& subMatrix6, Matrix<T>& subMatrix7, Matrix<T>& subMatrix8,
    Matrix<T>& subMatrix9, Matrix<T>& subMatrix10, Matrix<T>& subMatrix11, Matrix<T>& subMatrix12,
    Matrix<T>& subMatrix13, Matrix<T>& subMatrix14, Matrix<T>& subMatrix15, Matrix<T>& subMatrix16) const {
        std::size_t subSize = d;

        for (std::size_t i = 0; i < subSize; ++i) {
            for (std::size_t j = 0; j < subSize; ++j) {
                subMatrix1.set(j, i, M[i][j]);
                subMatrix2.set(j, i, M[i][j + subSize]);
                subMatrix3.set(j, i, M[i][j + 2*subSize]);
                subMatrix4.set(j, i, M[i][j + 3*subSize]);
                subMatrix5.set(j, i, M[i + subSize][j]);
                subMatrix6.set(j, i, M[i + subSize][j + subSize]);
                subMatrix7.set(j, i, M[i + subSize][j + 2*subSize]);
                subMatrix8.set(j, i, M[i + subSize][j + 3*subSize]);
                subMatrix9.set(j, i, M[i + 2*subSize][j]);
                subMatrix10.set(j, i, M[i + 2*subSize][j + subSize]);
                subMatrix11.set(j, i, M[i + 2*subSize][j + 2*subSize]);
                subMatrix12.set(j, i, M[i + 2*subSize][j + 3*subSize]);
                subMatrix13.set(j, i, M[i + 3*subSize][j]);
                subMatrix14.set(j, i, M[i + 3*subSize][j + subSize]);
                subMatrix15.set(j, i, M[i + 3*subSize][j + 2*subSize]);
                subMatrix16.set(j, i, M[i + 3*subSize][j + 3*subSize]);

            }
        }
    }

}; // end of class matrix

void random_pack_matrix_generator_new_encoding(std::size_t n, std::size_t m, Matrix<double> & A);
void random_pack_matrix_generator(std::size_t n, std::size_t m, Matrix<double> & A);
void random_block_matrix_generator(std::size_t n, std::size_t m, Matrix<double> & A);
void random_block_matrix_generator_512(std::size_t n, std::size_t m, Matrix<double> & A);
void random_matrix_generator(std::size_t n, std::size_t m, Matrix<double> & A);
void random_tae_matrix_generator(std::size_t n, std::size_t m, Matrix<double> & A);
void random_pack_triangle_matrix_generator(std::size_t n, std::size_t m, Matrix<double> & A);
void eye_matrix_generator(std::size_t n, std::size_t m, Matrix<double> & A);
void store_matrix(std::size_t n, std::size_t m, Matrix<double> & A);
void store_vector(std::vector<double>& vec);
void store_vector_new_encoding(std::vector<double>& vec);
void store_vector_triangle(std::vector<double>& vec);
seal::Ciphertext pack_change_cipher_to_matrix_A(
    seal::SEALContext context,
    seal::Ciphertext& cipher_matrix, 
    double scale, 
    seal::CKKSEncoder &encoder, 
    seal::Evaluator &evaluator, 
    seal::Decryptor &decryptor,
    seal::GaloisKeys &galois_keys, 
    seal::RelinKeys &relin_keys);
seal::Ciphertext pack_change_cipher_to_matrix_A_new_encoding(
    seal::SEALContext context,
    seal::Ciphertext& cipher_matrix, 
    double scale, 
    seal::CKKSEncoder &encoder, 
    seal::Evaluator &evaluator, 
    seal::Decryptor &decryptor,
    seal::GaloisKeys &galois_keys, 
    seal::RelinKeys &relin_keys,
    size_t d);
seal::Ciphertext pack_change_cipher_to_matrix_B(
    seal::SEALContext context,
    seal::Ciphertext& cipher_matrix, 
    double scale, 
    seal::CKKSEncoder &encoder, 
    seal::Evaluator &evaluator, 
    seal::Decryptor &decryptor,
    seal::GaloisKeys &galois_keys, 
    seal::RelinKeys &relin_keys);
seal::Ciphertext pack_change_cipher_to_matrix_B_new_encoding(
    seal::SEALContext context,
    seal::Ciphertext& cipher_matrix, 
    double scale, 
    seal::CKKSEncoder &encoder, 
    seal::Evaluator &evaluator, 
    seal::Decryptor &decryptor,
    seal::GaloisKeys &galois_keys, 
    seal::RelinKeys &relin_keys,
    size_t d);
seal::Ciphertext pack_encrypted_matrix_mul(
    seal::SEALContext context,
    seal::Ciphertext& cipher_A, 
    seal::Ciphertext& cipher_B, 
    double scale, 
    seal::CKKSEncoder &encoder, 
    seal::Evaluator &evaluator, 
    seal::Decryptor &decryptor,
    seal::GaloisKeys &galois_keys, 
    seal::RelinKeys &relin_keys);
seal::Ciphertext pack_encrypted_matrix_mul_new_encoding(
    seal::SEALContext context,
    seal::Ciphertext& cipher_A, 
    seal::Ciphertext& cipher_B, 
    double scale, 
    seal::CKKSEncoder &encoder, 
    seal::Evaluator &evaluator, 
    seal::Decryptor &decryptor,
    seal::GaloisKeys &galois_keys, 
    seal::RelinKeys &relin_keys,
    size_t d);
void decrypt_and_print(seal::Ciphertext& cipher, 
    seal::Decryptor& decryptor, 
    seal::CKKSEncoder& encoder); 
void decrypt_and_print_new_encoding(seal::Ciphertext& cipher, 
    seal::Decryptor& decryptor, 
    seal::CKKSEncoder& encoder,
    size_t d);
void decrypt_and_print_new_encoding_over(seal::Ciphertext& cipher, 
    seal::Decryptor& decryptor, 
    seal::CKKSEncoder& encoder,
    std::vector<double> result_I,
    size_t d);
seal::Ciphertext rotate_as_16(seal::Ciphertext& cipher_matrix,int d, double scale, 
    seal::CKKSEncoder &encoder, 
    seal::Evaluator &evaluator, 
    seal::Decryptor &decryptor,
    seal::GaloisKeys &galois_keys, 
    seal::RelinKeys &relin_keys);
std::vector<double> read_data_from_file(const std::string& file_path);
#endif 