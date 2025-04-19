#include <iostream>
#include <seal/seal.h>
#include "matrix.h"

using namespace std;
using namespace seal;

class cipher_matrix_jiang{
private:
    size_t n;//rows
    size_t m;//cols
    size_t d;
    seal::Ciphertext cipher_matrix;

public:
    cipher_matrix_jiang();
    cipher_matrix_jiang(seal::Ciphertext &cipher,size_t n,size_t m,size_t d);
    cipher_matrix_jiang(Matrix<double> data,double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor);
    // void eye_matrix(Matrix<double> data, size_t n, size_t m);
    void enc_matrix_cipher(vector<double> &data,double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor);
    void dec_matrix_cipher(Matrix<double> &destination,seal::CKKSEncoder &encoder,seal::Decryptor &decryptor);
    virtual ~cipher_matrix_jiang();
    // change to matrix A
    void change_cipher_to_matrix_A(double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator,
                                   seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys);
    // change to matrix B
    void change_cipher_to_matrix_B(double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator,
                                   seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys);
                                       //rotate matrix A
    void rotate_ctA(int i ,seal::Ciphertext& destination,seal::Evaluator &evaluator,
                    seal::GaloisKeys &galois_keys, seal::RelinKeys& relin_keys, seal::CKKSEncoder& encoder);
    //rotate matrix B
    void rotate_ctB(int i ,seal::Ciphertext& destination,seal::Evaluator &evaluator,
                    seal::GaloisKeys &galois_keys, seal::RelinKeys& relin_keys, seal::CKKSEncoder& encoder);
    size_t get_rows();
    size_t get_cols();
    size_t get_matrix_mul_size();   
    seal::Ciphertext get_cipher_matrix();   
    // 返回cipher_matrix的非常量引用
    seal::Ciphertext& get_cipher_matrix_ref() {
        return cipher_matrix;
    }
    void match_levels(SEALContext context,seal::Ciphertext& cipher1, seal::Ciphertext& cipher2, seal::Evaluator &evaluator);

};
void match_levels(SEALContext context, Ciphertext& cipher1, Ciphertext& cipher2, Evaluator& evaluator) {
    // 确保两个密文的尺度相同
    // cout<<"wssss"<<endl;
    // ios old_fmt(nullptr);
    // old_fmt.copyfmt(cout);
    // cout<<fixed<<setprecision(10);
    // cout<<"  + exact scale"<<cipher1.scale()<<endl;
    // cout<<"+    e"<<cipher2.scale()<<endl;
    // cout<<endl;
    // cout.copyfmt(old_fmt);
    parms_id_type last=cipher2.parms_id();
    cipher1.scale()=cipher2.scale();
    // evaluator.rescale_to_inplace(cipher1,last);
    // cipher2.scale()=pow(2.0,40);

    
    evaluator.mod_switch_to_inplace(cipher1,last);
    // evaluator.rescale_to_inplace(cipher1,last);

    // double scale1 = cipher1.scale();
    // double scale2 = cipher2.scale();
    // if (scale1 != scale2) {
    //     // 将两个密文的尺度调整到它们尺度中的最大值
    //     double target_scale = std::max(scale1, scale2);
    //     cout<<"begin"<<endl;
    //     // 这里需要一个循环，因为rescale_to_next_inplace会将尺度增加到下一个标准尺度，
    //     // 而不是直接设置到精确的尺度值
    //     while (cipher1.scale() < target_scale) {
    //         evaluator.rescale_to_next_inplace(cipher1);
    //     }
    //     cout<<"jieshu"<<endl;
    //     while (cipher2.scale() < target_scale) {
    //         evaluator.rescale_to_next_inplace(cipher2);
    //     }
    // }
    // // 获取两个密文的parms_id
    // parms_id_type parms_id1 = cipher1.parms_id();
    // parms_id_type parms_id2 = cipher2.parms_id();

    // // 比较两个密文的层级并使它们在同一层级
    // if (parms_id1 != parms_id2) {
    //     // 如果cipher1的层级高于cipher2，将cipher1降低到cipher2的层级
    //     if (context.get_context_data(parms_id1)->chain_index() > context.get_context_data(parms_id2)->chain_index()) {
    //         evaluator.mod_switch_to_inplace(cipher1, parms_id2);
    //     } else {
    //         // 如果cipher2的层级高于cipher1，将cipher2降低到cipher1的层级
    //         evaluator.mod_switch_to_inplace(cipher2, parms_id1);
    //     }
    // }
}

inline cipher_matrix_jiang::cipher_matrix_jiang()
{
    n=0;
    m=0;
}
inline cipher_matrix_jiang::cipher_matrix_jiang(seal::Ciphertext &cipher,size_t n,size_t m,size_t d)
{
    cipher_matrix=cipher;
    this->n=n;
    this->m=m;
    this->d=d;
}

// inline void cipher_matrix_jiang:: eye_matrix(Matrix<double> data, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor)
// {
//     this->n=data.get_rows();
//     this->m=data.get_cols();
//     d=sqrt(encoder.slot_count());
//     data.resize(d,d);
//     vector<double> flatten_data=data.flatten_matrix_to_rows_vector();
// }

inline cipher_matrix_jiang::cipher_matrix_jiang(Matrix<double> data, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor)
{
    this->n=data.get_rows();
    this->m=data.get_cols();
    d=sqrt(encoder.slot_count());
    data.resize(d,d);
    vector<double> flatten_data=data.flatten_matrix_to_rows_vector();
    enc_matrix_cipher(flatten_data,scale,encoder,encryptor);
}

inline void cipher_matrix_jiang::enc_matrix_cipher(vector<double> &data, double scale, seal::CKKSEncoder &encoder, seal::Encryptor &encryptor)
{
    if(data.size()>encoder.slot_count()){
        cerr << "  !!! the number of slot is not enough for the Matrix" << endl;
    }

    seal::Plaintext plain_tmp;
    encoder.encode(data,scale,plain_tmp);
    encryptor.encrypt(plain_tmp,cipher_matrix);
}

inline void cipher_matrix_jiang::dec_matrix_cipher(Matrix<double> &destination, seal::CKKSEncoder &encoder, seal::Decryptor &decryptor)
{
    //decrypte ciphertext
    Plaintext plain_tmp;
    vector<double> vec_tmp;
    decryptor.decrypt(cipher_matrix,plain_tmp);
    encoder.decode(plain_tmp,vec_tmp);

    destination.resize(d,d);
    for(size_t i=0;i<d;i++){
        for(size_t j=0;j<d;j++){
            destination.set(i,j,vec_tmp[i*d+j]);
        }
    }

    destination.resize(n,m);
}

inline cipher_matrix_jiang::~cipher_matrix_jiang()
{
    ;
}


inline void cipher_matrix_jiang::change_cipher_to_matrix_A(double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys)
{
    // cout<<"start change cipher to matrix A"<<endl;
    Matrix<double> u_sigma;
    u_sigma.generate_u_sigma(d, d);
    Ciphertext cipher_tmp;//save intermediate variables.
    vector<Ciphertext> cipher_rotate,cipher_result;//save rotate result vector
    Plaintext plain_tmp;
    vector<double> vec_tmp;
    Ciphertext destination;


    //start change matrixA
    // auto start = std::chrono::high_resolution_clock::now();
    //step1. determining the dimensions
    int rotate_size = 2 * d - 1;
    int rotate_outside_size = ceil(sqrt(rotate_size));
    int rotate_inside_size = ceil(double(rotate_size) / double(rotate_outside_size));
    // step2. pre-rotate ciphertext
    for (int i = 0; i < rotate_inside_size; i++) {
        evaluator.rotate_vector(cipher_matrix, -d+d*d+i+1, galois_keys, cipher_tmp);
        cipher_rotate.push_back(cipher_tmp);
    }

    //step3. starting step function 
    for (int i = 0; i < rotate_outside_size; i++) {
        vector<Ciphertext> cipher_vector;
        for (int j = 0; j < rotate_inside_size; j++) {
            vec_tmp = u_sigma.diag_vector(i * rotate_inside_size + j-d+1,encoder.slot_count());
            if (std::all_of(vec_tmp.begin(), vec_tmp.end(), [](double num) { return num == 0; }) == 1) {
                continue;
            }
            else {
                std::rotate(vec_tmp.rbegin(), vec_tmp.rbegin() + i * rotate_inside_size, vec_tmp.rend());
                encoder.encode(vec_tmp, cipher_matrix.scale(), plain_tmp);
                evaluator.mod_switch_to_inplace(plain_tmp, cipher_matrix.parms_id());
                evaluator.multiply_plain(cipher_rotate[j], plain_tmp, cipher_tmp);
                cipher_vector.push_back(cipher_tmp);
            }
        }
        // if(cipher_vector.size()==0){
        //     continue;
        // }
        evaluator.add_many(cipher_vector, cipher_tmp);
        evaluator.relinearize_inplace(cipher_tmp, relin_keys);
        evaluator.rescale_to_next_inplace(cipher_tmp);
        evaluator.rotate_vector_inplace(cipher_tmp, i * rotate_inside_size, galois_keys);
        cipher_result.push_back(cipher_tmp);
    }
    evaluator.add_many(cipher_result, destination);
    cipher_matrix=destination;
}

inline void cipher_matrix_jiang::change_cipher_to_matrix_B(double scale, seal::CKKSEncoder &encoder, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys)
{
    // cout<<"start change cipher to matrix B"<<endl;
    // generate u_tau matrix
    Matrix<double> u_tau;
    u_tau.generate_u_tau(d, d);
    Ciphertext cipher_tmp;//save intermediate variables
    vector<Ciphertext> cipher_rotate, cipher_result;//save rotate result vector
    Plaintext plain_tmp;
    vector<double> vec_tmp;
    Ciphertext destination;

    //start change matrixB
    //step1. determining the dimensions
    int rotate_size = d;
    int rotate_outside_size = ceil(sqrt(rotate_size));
    int rotate_inside_size = ceil(double(rotate_size) / double(rotate_outside_size));
    // step2. pre-rotate ciphertext 
    for (int i = 0; i < rotate_inside_size; i++) {
        evaluator.rotate_vector(cipher_matrix, i*d, galois_keys, cipher_tmp);
        cipher_rotate.push_back(cipher_tmp);
    }
    //step3. starting step function 
    for (int i = 0; i < rotate_outside_size; i++) {
        // cout<<i<<endl;
        vector<Ciphertext> cipher_vector;
        for (int j = 0; j < rotate_inside_size; j++) {
            if(size_t(i*rotate_inside_size+j)>=d){
                continue;
            }
            vec_tmp = u_tau.diag_vector((i*rotate_inside_size+j)*d,encoder.slot_count());
            // for(size_t k=0;k<512;k++){
            //     if(vec_tmp[k]==1){
            //         cout<<k<<"  ";
            //     }
            // }
            // cout<<endl;
            if (std::all_of(vec_tmp.begin(), vec_tmp.end(), [](double num) { return num == 0; }) == 1) {
                continue;
            }
            else {
                std::rotate(vec_tmp.rbegin(), vec_tmp.rbegin() + i * rotate_inside_size*d, vec_tmp.rend());
                encoder.encode(vec_tmp, cipher_matrix.scale(), plain_tmp);
                evaluator.mod_switch_to_inplace(plain_tmp, cipher_matrix.parms_id());
                evaluator.multiply_plain(cipher_rotate[j], plain_tmp, cipher_tmp);
                cipher_vector.push_back(cipher_tmp);
            }
        }
        // if(cipher_vector.size()==0){
        //     continue;
        // }
        evaluator.add_many(cipher_vector, cipher_tmp);
        evaluator.relinearize_inplace(cipher_tmp, relin_keys);
        evaluator.rescale_to_next_inplace(cipher_tmp);
        evaluator.rotate_vector_inplace(cipher_tmp, i * rotate_inside_size*d, galois_keys);
        cipher_result.push_back(cipher_tmp);
    }
    evaluator.add_many(cipher_result, destination);
    cipher_matrix=destination;
}

inline void cipher_matrix_jiang::rotate_ctA(int i, seal::Ciphertext &destination, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys, seal::CKKSEncoder &encoder)
{
    // step1.1.1 generate vector
    if (i != 0) {
        vector<double> rotate_right(d * d), rotate_left(d * d);
        seal::Plaintext right_plain, left_plain;
        seal::Ciphertext right_cipher, left_cipher, result_cipher;
        for (size_t j = 0; j < d; j++) {
            size_t k = 0;
            while (k < d) {
                if (k >= d-i) {
                    rotate_right[j * d + k] = 1;
                }
                else {
                    rotate_left[j * d + k] = 1;
                }
                k++;
            }
        }

        //step1.1.2 encode vector
        encoder.encode(rotate_right, cipher_matrix.scale(), right_plain);
        evaluator.mod_switch_to_inplace(right_plain, cipher_matrix.parms_id());
        encoder.encode(rotate_left, cipher_matrix.scale(), left_plain);
        evaluator.mod_switch_to_inplace(left_plain, cipher_matrix.parms_id());


        //step1.1.3 rotate and multiply
        evaluator.rotate_vector(cipher_matrix, d * d - d + i, galois_keys, right_cipher);
        evaluator.rotate_vector(cipher_matrix, i, galois_keys, left_cipher);
        evaluator.multiply_plain_inplace(right_cipher, right_plain);
        evaluator.multiply_plain_inplace(left_cipher, left_plain);
        evaluator.add(right_cipher, left_cipher, result_cipher);
        evaluator.relinearize_inplace(result_cipher, relin_keys);
        evaluator.rescale_to_next_inplace(result_cipher);
        destination = result_cipher;
    }
    else {
        vector<double> rotate_all(d * d, 1);
        Plaintext plain_tmp;
        Ciphertext result_cipher;
        encoder.encode(rotate_all, cipher_matrix.scale(), plain_tmp);
        evaluator.mod_switch_to_inplace(plain_tmp, cipher_matrix.parms_id());
        evaluator.multiply_plain(cipher_matrix, plain_tmp, result_cipher);
        evaluator.relinearize_inplace(result_cipher, relin_keys);
        evaluator.rescale_to_next_inplace(result_cipher);
        destination = result_cipher;
    }
}

inline void cipher_matrix_jiang::rotate_ctB(int i, seal::Ciphertext &destination, seal::Evaluator &evaluator, seal::GaloisKeys &galois_keys, seal::RelinKeys &relin_keys, seal::CKKSEncoder &encoder)
{
    evaluator.rotate_vector(cipher_matrix,i,galois_keys,destination);
}

inline size_t cipher_matrix_jiang::get_rows()
{
    return n;
}

inline size_t cipher_matrix_jiang::get_cols()
{
    return m;
}

inline size_t cipher_matrix_jiang::get_matrix_mul_size()
{
    return d;
}



inline seal::Ciphertext cipher_matrix_jiang::get_cipher_matrix()
{
    return cipher_matrix;
}

void encrypted_jiang_matrix_multiplication(seal::CKKSEncoder &encoder, // seal::Decryptor &decryptor,
                                          double scale, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys,
                                          seal::GaloisKeys &galois_keys, cipher_matrix_jiang &encrypted_A,
                                          cipher_matrix_jiang &encrypted_B, cipher_matrix_jiang &destination)
{
    if(encrypted_A.get_cols()!=encrypted_B.get_rows()){
        cerr << "Matrix dimensions are not equal to multiply.";
    }
    size_t n=encrypted_A.get_rows();
    size_t p=encrypted_B.get_cols();
    size_t d=encrypted_A.get_matrix_mul_size();

    vector<Ciphertext>  cipher_vector;
    Ciphertext cipher_tmpA,cipher_tmpB,cipher_tmp,encrypted_C;//save intermediate variables
    Plaintext plain_tmp;
    vector<double> vec_tmp;

    for (size_t i = 0; i < d; i++) {
        //step1.1 move ct.A(0)-->ct.A(i)
        encrypted_A.rotate_ctA(i, cipher_tmpA, evaluator, galois_keys, relin_keys, encoder);
        encrypted_B.rotate_ctB(d*i, cipher_tmpB, evaluator, galois_keys, relin_keys, encoder);
        evaluator.mod_switch_to_inplace(cipher_tmpB, cipher_tmpA.parms_id());
        evaluator.multiply(cipher_tmpA, cipher_tmpB, cipher_tmp);
        cipher_vector.push_back(cipher_tmp);
    }
    evaluator.add_many(cipher_vector, encrypted_C);
    evaluator.relinearize_inplace(encrypted_C, relin_keys);
    evaluator.rescale_to_next_inplace(encrypted_C);

    cipher_matrix_jiang ct_C(encrypted_C,n,p,d);
    destination=ct_C;
};
std::vector<Ciphertext> block_256_LeftMulRight_AddRight(SEALContext context,seal::CKKSEncoder &encoder, // seal::Decryptor &decryptor,
                                          double scale, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys,
                                          seal::GaloisKeys &gal_keys, cipher_matrix_jiang &cipher_A0,
                                          cipher_matrix_jiang &cipher_A1,cipher_matrix_jiang &cipher_A2,cipher_matrix_jiang &cipher_A3,
                                          cipher_matrix_jiang &cipher_B0,cipher_matrix_jiang &cipher_B1,
                                          cipher_matrix_jiang &cipher_B2,cipher_matrix_jiang &cipher_B3)
{
    std::vector<Ciphertext> destination;
    // seal::Ciphertext cipher_A0_ciphertext,cipher_A1_ciphertext,cipher_A2_ciphertext,cipher_A3_ciphertext;
    seal::Ciphertext cipher_B0_ciphertext,cipher_B1_ciphertext,cipher_B2_ciphertext,cipher_B3_ciphertext;
    // cipher_A0_ciphertext=cipher_A0.get_cipher_matrix();
    // cipher_A1_ciphertext=cipher_A1.get_cipher_matrix();
    // cipher_A2_ciphertext=cipher_A2.get_cipher_matrix();
    // cipher_A3_ciphertext=cipher_A3.get_cipher_matrix();
    cipher_B0_ciphertext=cipher_B0.get_cipher_matrix();
    cipher_B1_ciphertext=cipher_B1.get_cipher_matrix();
    cipher_B2_ciphertext=cipher_B2.get_cipher_matrix();
    cipher_B3_ciphertext=cipher_B3.get_cipher_matrix();
    cipher_matrix_jiang cipher_A00;
    cipher_matrix_jiang cipher_A11;
    cipher_matrix_jiang cipher_A22;
    cipher_matrix_jiang cipher_A33;
    cipher_matrix_jiang cipher_B00;
    cipher_matrix_jiang cipher_B22;
    cipher_matrix_jiang cipher_B11;
    cipher_matrix_jiang cipher_B33;
    cipher_B00=cipher_B0;
    cipher_B11=cipher_B1;
    cipher_B22=cipher_B2;
    cipher_B33=cipher_B3;
    cipher_A00=cipher_A0;
    cipher_A11=cipher_A1;
    cipher_A22=cipher_A2;
    cipher_A33=cipher_A3;

    cipher_A00.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B00.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A11.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B11.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A22.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B22.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A33.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B33.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    
    //compute A^2

    cipher_matrix_jiang cipher_A0B0;
    cipher_matrix_jiang cipher_A1B2;
    cipher_matrix_jiang cipher_A0B1;
    cipher_matrix_jiang cipher_A1B3;
    cipher_matrix_jiang cipher_A2B0;
    cipher_matrix_jiang cipher_A3B2;
    cipher_matrix_jiang cipher_A2B1;
    cipher_matrix_jiang cipher_A3B3;
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A00,cipher_B00,cipher_A0B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A11,cipher_B22,cipher_A1B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A00,cipher_B11,cipher_A0B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A11,cipher_B33,cipher_A1B3);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A22,cipher_B00,cipher_A2B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A33,cipher_B22,cipher_A3B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A22,cipher_B11,cipher_A2B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A33,cipher_B33,cipher_A3B3);
    
    seal::Ciphertext cipher_A0B0_ciphertext,cipher_A1B2_ciphertext,cipher_A0B1_ciphertext,cipher_A1B3_ciphertext,cipher_A2B0_ciphertext,cipher_A3B2_ciphertext,cipher_A2B1_ciphertext,cipher_A3B3_ciphertext;
    cipher_A0B0_ciphertext=cipher_A0B0.get_cipher_matrix();
    cipher_A1B2_ciphertext=cipher_A1B2.get_cipher_matrix();
    cipher_A0B1_ciphertext=cipher_A0B1.get_cipher_matrix();
    cipher_A1B3_ciphertext=cipher_A1B3.get_cipher_matrix();
    cipher_A2B0_ciphertext=cipher_A2B0.get_cipher_matrix();
    cipher_A3B2_ciphertext=cipher_A3B2.get_cipher_matrix();
    cipher_A2B1_ciphertext=cipher_A2B1.get_cipher_matrix();
    cipher_A3B3_ciphertext=cipher_A3B3.get_cipher_matrix();

    //compute A^2+A
    seal::Ciphertext cipher_temp_0,cipher_temp_1,cipher_temp_2,cipher_temp_3;
    evaluator.add(cipher_A0B0_ciphertext,cipher_A1B2_ciphertext,cipher_temp_0);
    evaluator.add(cipher_A0B1_ciphertext,cipher_A1B3_ciphertext,cipher_temp_1);
    evaluator.add(cipher_A2B0_ciphertext,cipher_A3B2_ciphertext,cipher_temp_2);
    evaluator.add(cipher_A2B1_ciphertext,cipher_A3B3_ciphertext,cipher_temp_3);
    seal::Ciphertext cipher_temp_0_A0,cipher_temp_1_A1,cipher_temp_2_A2,cipher_temp_3_A3;
    match_levels(context,cipher_B0_ciphertext,cipher_temp_0,evaluator);
    evaluator.add(cipher_temp_0,cipher_B0_ciphertext,cipher_temp_0_A0);

    match_levels(context,cipher_B1_ciphertext,cipher_temp_1,evaluator);
    evaluator.add(cipher_temp_1,cipher_B1_ciphertext,cipher_temp_1_A1);

    match_levels(context,cipher_B2_ciphertext,cipher_temp_2,evaluator);
    evaluator.add(cipher_temp_2,cipher_B2_ciphertext,cipher_temp_2_A2);

    match_levels(context,cipher_B3_ciphertext,cipher_temp_3,evaluator);
    evaluator.add(cipher_temp_3,cipher_B3_ciphertext,cipher_temp_3_A3);
    
    destination.push_back(cipher_temp_0_A0);
    destination.push_back(cipher_temp_1_A1);
    destination.push_back(cipher_temp_2_A2);
    destination.push_back(cipher_temp_3_A3);
    return destination;
}
std::vector<Ciphertext> block_256_LeftSquare(SEALContext context,seal::CKKSEncoder &encoder, //seal::Decryptor &decryptor,
                                          double scale, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys,
                                          seal::GaloisKeys &gal_keys, cipher_matrix_jiang &cipher_A0,
                                          cipher_matrix_jiang &cipher_A1,cipher_matrix_jiang &cipher_A2,cipher_matrix_jiang &cipher_A3,
                                          cipher_matrix_jiang &cipher_B0,cipher_matrix_jiang &cipher_B1,
                                          cipher_matrix_jiang &cipher_B2,cipher_matrix_jiang &cipher_B3)
{
    std::vector<Ciphertext> destination;
    cipher_matrix_jiang cipher_A00;
    cipher_matrix_jiang cipher_A11;
    cipher_matrix_jiang cipher_A22;
    cipher_matrix_jiang cipher_A33;
    cipher_matrix_jiang cipher_B00;
    cipher_matrix_jiang cipher_B22;
    cipher_matrix_jiang cipher_B11;
    cipher_matrix_jiang cipher_B33;
    cipher_B00=cipher_B0;
    cipher_B11=cipher_B1;
    cipher_B22=cipher_B2;
    cipher_B33=cipher_B3;
    cipher_A00=cipher_A0;
    cipher_A11=cipher_A1;
    cipher_A22=cipher_A2;
    cipher_A33=cipher_A3;

    cipher_A00.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B00.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A11.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B11.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A22.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B22.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A33.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B33.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    
    //compute A^2

    cipher_matrix_jiang cipher_A0B0;
    cipher_matrix_jiang cipher_A1B2;
    cipher_matrix_jiang cipher_A0B1;
    cipher_matrix_jiang cipher_A1B3;
    cipher_matrix_jiang cipher_A2B0;
    cipher_matrix_jiang cipher_A3B2;
    cipher_matrix_jiang cipher_A2B1;
    cipher_matrix_jiang cipher_A3B3;
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A00,cipher_B00,cipher_A0B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A11,cipher_B22,cipher_A1B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A00,cipher_B11,cipher_A0B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A11,cipher_B33,cipher_A1B3);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A22,cipher_B00,cipher_A2B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A33,cipher_B22,cipher_A3B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A22,cipher_B11,cipher_A2B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A33,cipher_B33,cipher_A3B3);
    // Matrix <double> w,s;
    // cout<<"A0BO"<<endl;
    // cipher_A0B0.dec_matrix_cipher(s,encoder,decryptor);
    
    // cipher_A0B1.dec_matrix_cipher(w,encoder,decryptor);
    // s.print();
    // cout<<"---------------------------------"<<endl;
    // cout<<"A0B1"<<endl;
    // w.print();
    
    seal::Ciphertext cipher_A0B0_ciphertext,cipher_A1B2_ciphertext,cipher_A0B1_ciphertext,cipher_A1B3_ciphertext,cipher_A2B0_ciphertext,cipher_A3B2_ciphertext,cipher_A2B1_ciphertext,cipher_A3B3_ciphertext;
    cipher_A0B0_ciphertext=cipher_A0B0.get_cipher_matrix();
    cipher_A1B2_ciphertext=cipher_A1B2.get_cipher_matrix();
    cipher_A0B1_ciphertext=cipher_A0B1.get_cipher_matrix();
    cipher_A1B3_ciphertext=cipher_A1B3.get_cipher_matrix();
    cipher_A2B0_ciphertext=cipher_A2B0.get_cipher_matrix();
    cipher_A3B2_ciphertext=cipher_A3B2.get_cipher_matrix();
    cipher_A2B1_ciphertext=cipher_A2B1.get_cipher_matrix();
    cipher_A3B3_ciphertext=cipher_A3B3.get_cipher_matrix();

    //compute A^2
    seal::Ciphertext cipher_temp_0,cipher_temp_1,cipher_temp_2,cipher_temp_3;
    evaluator.add(cipher_A0B0_ciphertext,cipher_A1B2_ciphertext,cipher_temp_0);
    evaluator.add(cipher_A0B1_ciphertext,cipher_A1B3_ciphertext,cipher_temp_1);
    evaluator.add(cipher_A2B0_ciphertext,cipher_A3B2_ciphertext,cipher_temp_2);
    evaluator.add(cipher_A2B1_ciphertext,cipher_A3B3_ciphertext,cipher_temp_3);
    
    destination.push_back(cipher_temp_0);
    destination.push_back(cipher_temp_1);
    destination.push_back(cipher_temp_2);
    destination.push_back(cipher_temp_3);
    return destination;
}
std::vector<Ciphertext> block_512_LeftSquare(SEALContext context,seal::CKKSEncoder &encoder, // seal::Decryptor &decryptor,
                                          double scale, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys,
                                          seal::GaloisKeys &gal_keys, cipher_matrix_jiang &cipher_A0,
                                          cipher_matrix_jiang &cipher_A1,cipher_matrix_jiang &cipher_A2,cipher_matrix_jiang &cipher_A3,
                                          cipher_matrix_jiang &cipher_A4,cipher_matrix_jiang &cipher_A5,cipher_matrix_jiang &cipher_A6,
                                          cipher_matrix_jiang &cipher_A7,cipher_matrix_jiang &cipher_A8,cipher_matrix_jiang &cipher_A9,
                                          cipher_matrix_jiang &cipher_A10,cipher_matrix_jiang &cipher_A11,cipher_matrix_jiang &cipher_A12,
                                          cipher_matrix_jiang &cipher_A13,cipher_matrix_jiang &cipher_A14,cipher_matrix_jiang &cipher_A15,
                                          cipher_matrix_jiang &cipher_B0,cipher_matrix_jiang &cipher_B1,
                                          cipher_matrix_jiang &cipher_B2,cipher_matrix_jiang &cipher_B3,cipher_matrix_jiang &cipher_B4,
                                          cipher_matrix_jiang &cipher_B5,cipher_matrix_jiang &cipher_B6,cipher_matrix_jiang &cipher_B7,
                                          cipher_matrix_jiang &cipher_B8,cipher_matrix_jiang &cipher_B9,cipher_matrix_jiang &cipher_B10,
                                          cipher_matrix_jiang &cipher_B11,cipher_matrix_jiang &cipher_B12,cipher_matrix_jiang &cipher_B13,
                                          cipher_matrix_jiang &cipher_B14,cipher_matrix_jiang &cipher_B15)
{
    std::vector<Ciphertext> destination;
    cipher_matrix_jiang cipher_A00;
    cipher_matrix_jiang cipher_A1_1;
    cipher_matrix_jiang cipher_A22;
    cipher_matrix_jiang cipher_A33;
    cipher_matrix_jiang cipher_A44;
    cipher_matrix_jiang cipher_A55;
    cipher_matrix_jiang cipher_A66;
    cipher_matrix_jiang cipher_A77;
    cipher_matrix_jiang cipher_A88;
    cipher_matrix_jiang cipher_A99;
    cipher_matrix_jiang cipher_A1010;
    cipher_matrix_jiang cipher_A1111;
    cipher_matrix_jiang cipher_A1212;
    cipher_matrix_jiang cipher_A1313;
    cipher_matrix_jiang cipher_A1414;
    cipher_matrix_jiang cipher_A1515;

    cipher_matrix_jiang cipher_B00;
    cipher_matrix_jiang cipher_B1_1;
    cipher_matrix_jiang cipher_B22;
    cipher_matrix_jiang cipher_B33;
    cipher_matrix_jiang cipher_B44;
    cipher_matrix_jiang cipher_B55;
    cipher_matrix_jiang cipher_B66;
    cipher_matrix_jiang cipher_B77;
    cipher_matrix_jiang cipher_B88;
    cipher_matrix_jiang cipher_B99;
    cipher_matrix_jiang cipher_B1010;
    cipher_matrix_jiang cipher_B1111;
    cipher_matrix_jiang cipher_B1212;
    cipher_matrix_jiang cipher_B1313;
    cipher_matrix_jiang cipher_B1414;
    cipher_matrix_jiang cipher_B1515;
    
    cipher_A00=cipher_A0;
    cipher_A1_1=cipher_A1;
    cipher_A22=cipher_A2;
    cipher_A33=cipher_A3;
    cipher_A44=cipher_A4;
    cipher_A55=cipher_A5;
    cipher_A66=cipher_A6;
    cipher_A77=cipher_A7;
    cipher_A88=cipher_A8;
    cipher_A99=cipher_A9;
    cipher_A1010=cipher_A10;
    cipher_A1111=cipher_A11;
    cipher_A1212=cipher_A12;
    cipher_A1313=cipher_A13;
    cipher_A1414=cipher_A14;
    cipher_A1515=cipher_A15;
    
    cipher_B00=cipher_B0;
    cipher_B1_1=cipher_B1;
    cipher_B22=cipher_B2;
    cipher_B33=cipher_B3;
    cipher_B44=cipher_B4;
    cipher_B55=cipher_B5;
    cipher_B66=cipher_B6;
    cipher_B77=cipher_B7;
    cipher_B88=cipher_B8;
    cipher_B99=cipher_B9;
    cipher_B1010=cipher_B10;
    cipher_B1111=cipher_B11;
    cipher_B1212=cipher_B12;
    cipher_B1313=cipher_B13;
    cipher_B1414=cipher_B14;
    cipher_B1515=cipher_B15;

    cipher_A00.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B00.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A1_1.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B1_1.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A22.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B22.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A33.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B33.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A44.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B44.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A55.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B55.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A66.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B66.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A77.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B77.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A88.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B88.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A99.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B99.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A1010.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B1010.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A1111.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B1111.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A1212.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B1212.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A1313.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B1313.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A1414.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B1414.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A1515.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B1515.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    
    cipher_matrix_jiang cipher_A0B0,cipher_A1B4,cipher_A2B8,cipher_A3B12;
    cipher_matrix_jiang cipher_A0B1,cipher_A1B5,cipher_A2B9,cipher_A3B13;
    cipher_matrix_jiang cipher_A0B2,cipher_A1B6,cipher_A2B10,cipher_A3B14;
    cipher_matrix_jiang cipher_A0B3,cipher_A1B7,cipher_A2B11,cipher_A3B15;

    cipher_matrix_jiang cipher_A4B0,cipher_A5B4,cipher_A6B8,cipher_A7B12;
    cipher_matrix_jiang cipher_A4B1,cipher_A5B5,cipher_A6B9,cipher_A7B13;
    cipher_matrix_jiang cipher_A4B2,cipher_A5B6,cipher_A6B10,cipher_A7B14;
    cipher_matrix_jiang cipher_A4B3,cipher_A5B7,cipher_A6B11,cipher_A7B15;

    cipher_matrix_jiang cipher_A8B0,cipher_A9B4,cipher_A10B8,cipher_A11B12;
    cipher_matrix_jiang cipher_A8B1,cipher_A9B5,cipher_A10B9,cipher_A11B13;
    cipher_matrix_jiang cipher_A8B2,cipher_A9B6,cipher_A10B10,cipher_A11B14;
    cipher_matrix_jiang cipher_A8B3,cipher_A9B7,cipher_A10B11,cipher_A11B15;

    cipher_matrix_jiang cipher_A12B0,cipher_A13B4,cipher_A14B8,cipher_A15B12;
    cipher_matrix_jiang cipher_A12B1,cipher_A13B5,cipher_A14B9,cipher_A15B13;
    cipher_matrix_jiang cipher_A12B2,cipher_A13B6,cipher_A14B10,cipher_A15B14;
    cipher_matrix_jiang cipher_A12B3,cipher_A13B7,cipher_A14B11,cipher_A15B15;

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A00,cipher_B00,cipher_A0B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A00,cipher_B1_1,cipher_A0B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A00,cipher_B22,cipher_A0B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A00,cipher_B33,cipher_A0B3);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A44,cipher_B00,cipher_A4B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A44,cipher_B1_1,cipher_A4B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A44,cipher_B22,cipher_A4B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A44,cipher_B33,cipher_A4B3);
        
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A88,cipher_B00,cipher_A8B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A88,cipher_B1_1,cipher_A8B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A88,cipher_B22,cipher_A8B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A88,cipher_B33,cipher_A8B3);
        
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1212,cipher_B00,cipher_A12B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1212,cipher_B1_1,cipher_A12B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1212,cipher_B22,cipher_A12B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1212,cipher_B33,cipher_A12B3);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1_1,cipher_B44,cipher_A1B4);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1_1,cipher_B55,cipher_A1B5);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1_1,cipher_B66,cipher_A1B6);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1_1,cipher_B77,cipher_A1B7);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A55,cipher_B44,cipher_A5B4);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A55,cipher_B55,cipher_A5B5);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A55,cipher_B66,cipher_A5B6);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A55,cipher_B77,cipher_A5B7);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A99,cipher_B44,cipher_A9B4);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A99,cipher_B55,cipher_A9B5);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A99,cipher_B66,cipher_A9B6);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A99,cipher_B77,cipher_A9B7);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1313,cipher_B44,cipher_A13B4);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1313,cipher_B55,cipher_A13B5);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1313,cipher_B66,cipher_A13B6);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1313,cipher_B77,cipher_A13B7);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A22,cipher_B88,cipher_A2B8);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A22,cipher_B99,cipher_A2B9);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A22,cipher_B1010,cipher_A2B10);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A22,cipher_B1111,cipher_A2B11);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A66,cipher_B88,cipher_A6B8);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A66,cipher_B99,cipher_A6B9);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A66,cipher_B1010,cipher_A6B10);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A66,cipher_B1111,cipher_A6B11);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1010,cipher_B88,cipher_A10B8);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1010,cipher_B99,cipher_A10B9);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1010,cipher_B1010,cipher_A10B10);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1010,cipher_B1111,cipher_A10B11);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1414,cipher_B88,cipher_A14B8);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1414,cipher_B99,cipher_A14B9);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1414,cipher_B1010,cipher_A14B10);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1414,cipher_B1111,cipher_A14B11);
       
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A33,cipher_B1212,cipher_A3B12);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A33,cipher_B1313,cipher_A3B13);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A33,cipher_B1414,cipher_A3B14);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A33,cipher_B1515,cipher_A3B15);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A77,cipher_B1212,cipher_A7B12);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A77,cipher_B1313,cipher_A7B13);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A77,cipher_B1414,cipher_A7B14);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A77,cipher_B1515,cipher_A7B15);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1111,cipher_B1212,cipher_A11B12);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1111,cipher_B1313,cipher_A11B13);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1111,cipher_B1414,cipher_A11B14);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1111,cipher_B1515,cipher_A11B15);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1515,cipher_B1212,cipher_A15B12);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1515,cipher_B1313,cipher_A15B13);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1515,cipher_B1414,cipher_A15B14);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1515,cipher_B1515,cipher_A15B15);
   
    seal::Ciphertext cipher_A0B0_ciphertext,cipher_A1B4_ciphertext,cipher_A2B8_ciphertext,cipher_A3B12_ciphertext;
    seal::Ciphertext cipher_A0B1_ciphertext,cipher_A1B5_ciphertext,cipher_A2B9_ciphertext,cipher_A3B13_ciphertext;
    seal::Ciphertext cipher_A0B2_ciphertext,cipher_A1B6_ciphertext,cipher_A2B10_ciphertext,cipher_A3B14_ciphertext;
    seal::Ciphertext cipher_A0B3_ciphertext,cipher_A1B7_ciphertext,cipher_A2B11_ciphertext,cipher_A3B15_ciphertext;
    seal::Ciphertext cipher_A4B0_ciphertext,cipher_A5B4_ciphertext,cipher_A6B8_ciphertext,cipher_A7B12_ciphertext;
    seal::Ciphertext cipher_A4B1_ciphertext,cipher_A5B5_ciphertext,cipher_A6B9_ciphertext,cipher_A7B13_ciphertext;
    seal::Ciphertext cipher_A4B2_ciphertext,cipher_A5B6_ciphertext,cipher_A6B10_ciphertext,cipher_A7B14_ciphertext;
    seal::Ciphertext cipher_A4B3_ciphertext,cipher_A5B7_ciphertext,cipher_A6B11_ciphertext,cipher_A7B15_ciphertext;

    seal::Ciphertext cipher_A8B0_ciphertext,cipher_A9B4_ciphertext,cipher_A10B8_ciphertext,cipher_A11B12_ciphertext;
    seal::Ciphertext cipher_A8B1_ciphertext,cipher_A9B5_ciphertext,cipher_A10B9_ciphertext,cipher_A11B13_ciphertext;
    seal::Ciphertext cipher_A8B2_ciphertext,cipher_A9B6_ciphertext,cipher_A10B10_ciphertext,cipher_A11B14_ciphertext;
    seal::Ciphertext cipher_A8B3_ciphertext,cipher_A9B7_ciphertext,cipher_A10B11_ciphertext,cipher_A11B15_ciphertext;

    seal::Ciphertext cipher_A12B0_ciphertext,cipher_A13B4_ciphertext,cipher_A14B8_ciphertext,cipher_A15B12_ciphertext;
    seal::Ciphertext cipher_A12B1_ciphertext,cipher_A13B5_ciphertext,cipher_A14B9_ciphertext,cipher_A15B13_ciphertext;
    seal::Ciphertext cipher_A12B2_ciphertext,cipher_A13B6_ciphertext,cipher_A14B10_ciphertext,cipher_A15B14_ciphertext;
    seal::Ciphertext cipher_A12B3_ciphertext,cipher_A13B7_ciphertext,cipher_A14B11_ciphertext,cipher_A15B15_ciphertext;
    
    seal::Ciphertext cipher_temp_0,cipher_temp_1,cipher_temp_2,cipher_temp_3;
    seal::Ciphertext cipher_temp_4,cipher_temp_5,cipher_temp_6,cipher_temp_7;
    seal::Ciphertext cipher_temp_8,cipher_temp_9,cipher_temp_10,cipher_temp_11;
    seal::Ciphertext cipher_temp_12,cipher_temp_13,cipher_temp_14,cipher_temp_15;

    cipher_A0B0_ciphertext=cipher_A0B0.get_cipher_matrix();
    cipher_A0B1_ciphertext=cipher_A0B1.get_cipher_matrix();
    cipher_A0B2_ciphertext=cipher_A0B2.get_cipher_matrix();
    cipher_A0B3_ciphertext=cipher_A0B3.get_cipher_matrix();
    cipher_A4B0_ciphertext=cipher_A4B0.get_cipher_matrix();
    cipher_A4B1_ciphertext=cipher_A4B1.get_cipher_matrix();
    cipher_A4B2_ciphertext=cipher_A4B2.get_cipher_matrix();
    cipher_A4B3_ciphertext=cipher_A4B3.get_cipher_matrix();
    cipher_A8B0_ciphertext=cipher_A8B0.get_cipher_matrix();
    cipher_A8B1_ciphertext=cipher_A8B1.get_cipher_matrix();
    cipher_A8B2_ciphertext=cipher_A8B2.get_cipher_matrix();
    cipher_A8B3_ciphertext=cipher_A8B3.get_cipher_matrix();
    cipher_A12B0_ciphertext=cipher_A12B0.get_cipher_matrix();
    cipher_A12B1_ciphertext=cipher_A12B1.get_cipher_matrix();
    cipher_A12B2_ciphertext=cipher_A12B2.get_cipher_matrix();
    cipher_A12B3_ciphertext=cipher_A12B3.get_cipher_matrix();

    cipher_A1B4_ciphertext=cipher_A1B4.get_cipher_matrix();
    cipher_A1B5_ciphertext=cipher_A1B5.get_cipher_matrix();
    cipher_A1B6_ciphertext=cipher_A1B6.get_cipher_matrix();
    cipher_A1B7_ciphertext=cipher_A1B7.get_cipher_matrix();
    cipher_A5B4_ciphertext=cipher_A5B4.get_cipher_matrix();
    cipher_A5B5_ciphertext=cipher_A5B5.get_cipher_matrix();
    cipher_A5B6_ciphertext=cipher_A5B6.get_cipher_matrix();
    cipher_A5B7_ciphertext=cipher_A5B7.get_cipher_matrix();
    cipher_A9B4_ciphertext=cipher_A9B4.get_cipher_matrix();
    cipher_A9B5_ciphertext=cipher_A9B5.get_cipher_matrix();
    cipher_A9B6_ciphertext=cipher_A9B6.get_cipher_matrix();
    cipher_A9B7_ciphertext=cipher_A9B7.get_cipher_matrix();
    cipher_A13B4_ciphertext=cipher_A13B4.get_cipher_matrix();
    cipher_A13B5_ciphertext=cipher_A13B5.get_cipher_matrix();
    cipher_A13B6_ciphertext=cipher_A13B6.get_cipher_matrix();
    cipher_A13B7_ciphertext=cipher_A13B7.get_cipher_matrix();

    cipher_A2B8_ciphertext=cipher_A2B8.get_cipher_matrix();
    cipher_A2B9_ciphertext=cipher_A2B9.get_cipher_matrix();
    cipher_A2B10_ciphertext=cipher_A2B10.get_cipher_matrix();
    cipher_A2B11_ciphertext=cipher_A2B11.get_cipher_matrix();
    cipher_A6B8_ciphertext=cipher_A6B8.get_cipher_matrix();
    cipher_A6B9_ciphertext=cipher_A6B9.get_cipher_matrix();
    cipher_A6B10_ciphertext=cipher_A6B10.get_cipher_matrix();
    cipher_A6B11_ciphertext=cipher_A6B11.get_cipher_matrix();
    cipher_A10B8_ciphertext=cipher_A10B8.get_cipher_matrix();
    cipher_A10B9_ciphertext=cipher_A10B9.get_cipher_matrix();
    cipher_A10B10_ciphertext=cipher_A10B10.get_cipher_matrix();
    cipher_A10B11_ciphertext=cipher_A10B11.get_cipher_matrix();
    cipher_A14B8_ciphertext=cipher_A14B8.get_cipher_matrix();
    cipher_A14B9_ciphertext=cipher_A14B9.get_cipher_matrix();
    cipher_A14B10_ciphertext=cipher_A14B10.get_cipher_matrix();
    cipher_A14B11_ciphertext=cipher_A14B11.get_cipher_matrix();

    cipher_A3B12_ciphertext=cipher_A3B12.get_cipher_matrix();
    cipher_A3B13_ciphertext=cipher_A3B13.get_cipher_matrix();
    cipher_A3B14_ciphertext=cipher_A3B14.get_cipher_matrix();
    cipher_A3B15_ciphertext=cipher_A3B15.get_cipher_matrix();
    cipher_A7B12_ciphertext=cipher_A7B12.get_cipher_matrix();
    cipher_A7B13_ciphertext=cipher_A7B13.get_cipher_matrix();
    cipher_A7B14_ciphertext=cipher_A7B14.get_cipher_matrix();
    cipher_A7B15_ciphertext=cipher_A7B15.get_cipher_matrix();
    cipher_A11B12_ciphertext=cipher_A11B12.get_cipher_matrix();
    cipher_A11B13_ciphertext=cipher_A11B13.get_cipher_matrix();
    cipher_A11B14_ciphertext=cipher_A11B14.get_cipher_matrix();
    cipher_A11B15_ciphertext=cipher_A11B15.get_cipher_matrix();
    cipher_A15B12_ciphertext=cipher_A15B12.get_cipher_matrix();
    cipher_A15B13_ciphertext=cipher_A15B13.get_cipher_matrix();
    cipher_A15B14_ciphertext=cipher_A15B14.get_cipher_matrix();
    cipher_A15B15_ciphertext=cipher_A15B15.get_cipher_matrix();

    evaluator.add(cipher_A0B0_ciphertext,cipher_A1B4_ciphertext,cipher_temp_0);
    evaluator.add_inplace(cipher_temp_0,cipher_A2B8_ciphertext);
    evaluator.add_inplace(cipher_temp_0,cipher_A3B12_ciphertext);
    evaluator.add(cipher_A0B1_ciphertext,cipher_A1B5_ciphertext,cipher_temp_1);
    evaluator.add_inplace(cipher_temp_1,cipher_A2B9_ciphertext);
    evaluator.add_inplace(cipher_temp_1,cipher_A3B13_ciphertext);
    evaluator.add(cipher_A0B2_ciphertext,cipher_A1B6_ciphertext,cipher_temp_2);
    evaluator.add_inplace(cipher_temp_2,cipher_A2B10_ciphertext);
    evaluator.add_inplace(cipher_temp_2,cipher_A3B14_ciphertext);
    evaluator.add(cipher_A0B3_ciphertext,cipher_A1B7_ciphertext,cipher_temp_3);
    evaluator.add_inplace(cipher_temp_3,cipher_A2B11_ciphertext);
    evaluator.add_inplace(cipher_temp_3,cipher_A3B15_ciphertext);

    evaluator.add(cipher_A4B0_ciphertext,cipher_A5B4_ciphertext,cipher_temp_4);
    evaluator.add_inplace(cipher_temp_4,cipher_A6B8_ciphertext);
    evaluator.add_inplace(cipher_temp_4,cipher_A7B12_ciphertext);
    evaluator.add(cipher_A4B1_ciphertext,cipher_A5B5_ciphertext,cipher_temp_5);
    evaluator.add_inplace(cipher_temp_5,cipher_A6B9_ciphertext);
    evaluator.add_inplace(cipher_temp_5,cipher_A7B13_ciphertext);
    evaluator.add(cipher_A4B2_ciphertext,cipher_A5B6_ciphertext,cipher_temp_6);
    evaluator.add_inplace(cipher_temp_6,cipher_A6B10_ciphertext);
    evaluator.add_inplace(cipher_temp_6,cipher_A7B14_ciphertext);
    evaluator.add(cipher_A4B3_ciphertext,cipher_A5B7_ciphertext,cipher_temp_7);
    evaluator.add_inplace(cipher_temp_7,cipher_A6B11_ciphertext);
    evaluator.add_inplace(cipher_temp_7,cipher_A7B15_ciphertext);

    evaluator.add(cipher_A8B0_ciphertext,cipher_A9B4_ciphertext,cipher_temp_8);
    evaluator.add_inplace(cipher_temp_8,cipher_A10B8_ciphertext);
    evaluator.add_inplace(cipher_temp_8,cipher_A11B12_ciphertext);
    evaluator.add(cipher_A8B1_ciphertext,cipher_A9B5_ciphertext,cipher_temp_9);
    evaluator.add_inplace(cipher_temp_9,cipher_A10B9_ciphertext);
    evaluator.add_inplace(cipher_temp_9,cipher_A11B13_ciphertext);
    evaluator.add(cipher_A8B2_ciphertext,cipher_A9B6_ciphertext,cipher_temp_10);
    evaluator.add_inplace(cipher_temp_10,cipher_A10B10_ciphertext);
    evaluator.add_inplace(cipher_temp_10,cipher_A11B14_ciphertext);
    evaluator.add(cipher_A8B3_ciphertext,cipher_A9B7_ciphertext,cipher_temp_11);
    evaluator.add_inplace(cipher_temp_11,cipher_A10B11_ciphertext);
    evaluator.add_inplace(cipher_temp_11,cipher_A11B15_ciphertext);

    evaluator.add(cipher_A12B0_ciphertext,cipher_A13B4_ciphertext,cipher_temp_12);
    evaluator.add_inplace(cipher_temp_12,cipher_A14B8_ciphertext);
    evaluator.add_inplace(cipher_temp_12,cipher_A15B12_ciphertext);
    evaluator.add(cipher_A12B1_ciphertext,cipher_A13B5_ciphertext,cipher_temp_13);
    evaluator.add_inplace(cipher_temp_13,cipher_A14B9_ciphertext);
    evaluator.add_inplace(cipher_temp_13,cipher_A15B13_ciphertext);
    evaluator.add(cipher_A12B2_ciphertext,cipher_A13B6_ciphertext,cipher_temp_14);
    evaluator.add_inplace(cipher_temp_14,cipher_A14B10_ciphertext);
    evaluator.add_inplace(cipher_temp_14,cipher_A15B14_ciphertext);
    evaluator.add(cipher_A12B3_ciphertext,cipher_A13B7_ciphertext,cipher_temp_15);
    evaluator.add_inplace(cipher_temp_15,cipher_A14B11_ciphertext);
    evaluator.add_inplace(cipher_temp_15,cipher_A15B15_ciphertext);
    destination.push_back(cipher_temp_0);
    destination.push_back(cipher_temp_1);
    destination.push_back(cipher_temp_2);
    destination.push_back(cipher_temp_3);
    destination.push_back(cipher_temp_4);
    destination.push_back(cipher_temp_5);
    destination.push_back(cipher_temp_6);
    destination.push_back(cipher_temp_7);
    destination.push_back(cipher_temp_8);
    destination.push_back(cipher_temp_9);
    destination.push_back(cipher_temp_10);
    destination.push_back(cipher_temp_11);
    destination.push_back(cipher_temp_12);
    destination.push_back(cipher_temp_13);
    destination.push_back(cipher_temp_14);
    destination.push_back(cipher_temp_15);
    return destination;
}
std::vector<Ciphertext> block_512_LeftMulRight_AddRight(SEALContext context,seal::CKKSEncoder &encoder,
                                          double scale, seal::Evaluator &evaluator, seal::RelinKeys &relin_keys,
                                          seal::GaloisKeys &gal_keys, cipher_matrix_jiang &cipher_A0,
                                          cipher_matrix_jiang &cipher_A1,cipher_matrix_jiang &cipher_A2,cipher_matrix_jiang &cipher_A3,
                                          cipher_matrix_jiang &cipher_A4,cipher_matrix_jiang &cipher_A5,cipher_matrix_jiang &cipher_A6,
                                          cipher_matrix_jiang &cipher_A7,cipher_matrix_jiang &cipher_A8,cipher_matrix_jiang &cipher_A9,
                                          cipher_matrix_jiang &cipher_A10,cipher_matrix_jiang &cipher_A11,cipher_matrix_jiang &cipher_A12,
                                          cipher_matrix_jiang &cipher_A13,cipher_matrix_jiang &cipher_A14,cipher_matrix_jiang &cipher_A15,
                                          cipher_matrix_jiang &cipher_B0,cipher_matrix_jiang &cipher_B1,
                                          cipher_matrix_jiang &cipher_B2,cipher_matrix_jiang &cipher_B3,cipher_matrix_jiang &cipher_B4,
                                          cipher_matrix_jiang &cipher_B5,cipher_matrix_jiang &cipher_B6,cipher_matrix_jiang &cipher_B7,
                                          cipher_matrix_jiang &cipher_B8,cipher_matrix_jiang &cipher_B9,cipher_matrix_jiang &cipher_B10,
                                          cipher_matrix_jiang &cipher_B11,cipher_matrix_jiang &cipher_B12,cipher_matrix_jiang &cipher_B13,
                                          cipher_matrix_jiang &cipher_B14,cipher_matrix_jiang &cipher_B15)
{
    std::vector<Ciphertext> destination;
    seal::Ciphertext cipher_A0_ciphertext,cipher_A1_ciphertext,cipher_A2_ciphertext,cipher_A3_ciphertext;
    seal::Ciphertext cipher_A4_ciphertext,cipher_A5_ciphertext,cipher_A6_ciphertext,cipher_A7_ciphertext;
    seal::Ciphertext cipher_A8_ciphertext,cipher_A9_ciphertext,cipher_A10_ciphertext,cipher_A11_ciphertext;
    seal::Ciphertext cipher_A12_ciphertext,cipher_A13_ciphertext,cipher_A14_ciphertext,cipher_A15_ciphertext;
 

    cipher_A0_ciphertext=cipher_B0.get_cipher_matrix();
    cipher_A1_ciphertext=cipher_B1.get_cipher_matrix();
    cipher_A2_ciphertext=cipher_B2.get_cipher_matrix();
    cipher_A3_ciphertext=cipher_B3.get_cipher_matrix();
    cipher_A4_ciphertext=cipher_B4.get_cipher_matrix();
    cipher_A5_ciphertext=cipher_B5.get_cipher_matrix();
    cipher_A6_ciphertext=cipher_B6.get_cipher_matrix();
    cipher_A7_ciphertext=cipher_B7.get_cipher_matrix();
    cipher_A8_ciphertext=cipher_B8.get_cipher_matrix();
    cipher_A9_ciphertext=cipher_B9.get_cipher_matrix();
    cipher_A10_ciphertext=cipher_B10.get_cipher_matrix();
    cipher_A11_ciphertext=cipher_B11.get_cipher_matrix();
    cipher_A12_ciphertext=cipher_B12.get_cipher_matrix();
    cipher_A13_ciphertext=cipher_B13.get_cipher_matrix();
    cipher_A14_ciphertext=cipher_B14.get_cipher_matrix();
    cipher_A15_ciphertext=cipher_B15.get_cipher_matrix();

    cipher_matrix_jiang cipher_A00;
    cipher_matrix_jiang cipher_A1_1;
    cipher_matrix_jiang cipher_A22;
    cipher_matrix_jiang cipher_A33;
    cipher_matrix_jiang cipher_A44;
    cipher_matrix_jiang cipher_A55;
    cipher_matrix_jiang cipher_A66;
    cipher_matrix_jiang cipher_A77;
    cipher_matrix_jiang cipher_A88;
    cipher_matrix_jiang cipher_A99;
    cipher_matrix_jiang cipher_A1010;
    cipher_matrix_jiang cipher_A1111;
    cipher_matrix_jiang cipher_A1212;
    cipher_matrix_jiang cipher_A1313;
    cipher_matrix_jiang cipher_A1414;
    cipher_matrix_jiang cipher_A1515;

    cipher_matrix_jiang cipher_B00;
    cipher_matrix_jiang cipher_B1_1;
    cipher_matrix_jiang cipher_B22;
    cipher_matrix_jiang cipher_B33;
    cipher_matrix_jiang cipher_B44;
    cipher_matrix_jiang cipher_B55;
    cipher_matrix_jiang cipher_B66;
    cipher_matrix_jiang cipher_B77;
    cipher_matrix_jiang cipher_B88;
    cipher_matrix_jiang cipher_B99;
    cipher_matrix_jiang cipher_B1010;
    cipher_matrix_jiang cipher_B1111;
    cipher_matrix_jiang cipher_B1212;
    cipher_matrix_jiang cipher_B1313;
    cipher_matrix_jiang cipher_B1414;
    cipher_matrix_jiang cipher_B1515;
    
    cipher_A00=cipher_A0;
    cipher_A1_1=cipher_A1;
    cipher_A22=cipher_A2;
    cipher_A33=cipher_A3;
    cipher_A44=cipher_A4;
    cipher_A55=cipher_A5;
    cipher_A66=cipher_A6;
    cipher_A77=cipher_A7;
    cipher_A88=cipher_A8;
    cipher_A99=cipher_A9;
    cipher_A1010=cipher_A10;
    cipher_A1111=cipher_A11;
    cipher_A1212=cipher_A12;
    cipher_A1313=cipher_A13;
    cipher_A1414=cipher_A14;
    cipher_A1515=cipher_A15;
    
    cipher_B00=cipher_B0;
    cipher_B1_1=cipher_B1;
    cipher_B22=cipher_B2;
    cipher_B33=cipher_B3;
    cipher_B44=cipher_B4;
    cipher_B55=cipher_B5;
    cipher_B66=cipher_B6;
    cipher_B77=cipher_B7;
    cipher_B88=cipher_B8;
    cipher_B99=cipher_B9;
    cipher_B1010=cipher_B10;
    cipher_B1111=cipher_B11;
    cipher_B1212=cipher_B12;
    cipher_B1313=cipher_B13;
    cipher_B1414=cipher_B14;
    cipher_B1515=cipher_B15;

    cipher_A00.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B00.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A1_1.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B1_1.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A22.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B22.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A33.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B33.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A44.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B44.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A55.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B55.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A66.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B66.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A77.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B77.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A88.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B88.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A99.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B99.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A1010.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B1010.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A1111.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B1111.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A1212.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B1212.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A1313.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B1313.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A1414.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B1414.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A1515.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B1515.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    
    cipher_matrix_jiang cipher_A0B0,cipher_A1B4,cipher_A2B8,cipher_A3B12;
    cipher_matrix_jiang cipher_A0B1,cipher_A1B5,cipher_A2B9,cipher_A3B13;
    cipher_matrix_jiang cipher_A0B2,cipher_A1B6,cipher_A2B10,cipher_A3B14;
    cipher_matrix_jiang cipher_A0B3,cipher_A1B7,cipher_A2B11,cipher_A3B15;

    cipher_matrix_jiang cipher_A4B0,cipher_A5B4,cipher_A6B8,cipher_A7B12;
    cipher_matrix_jiang cipher_A4B1,cipher_A5B5,cipher_A6B9,cipher_A7B13;
    cipher_matrix_jiang cipher_A4B2,cipher_A5B6,cipher_A6B10,cipher_A7B14;
    cipher_matrix_jiang cipher_A4B3,cipher_A5B7,cipher_A6B11,cipher_A7B15;

    cipher_matrix_jiang cipher_A8B0,cipher_A9B4,cipher_A10B8,cipher_A11B12;
    cipher_matrix_jiang cipher_A8B1,cipher_A9B5,cipher_A10B9,cipher_A11B13;
    cipher_matrix_jiang cipher_A8B2,cipher_A9B6,cipher_A10B10,cipher_A11B14;
    cipher_matrix_jiang cipher_A8B3,cipher_A9B7,cipher_A10B11,cipher_A11B15;

    cipher_matrix_jiang cipher_A12B0,cipher_A13B4,cipher_A14B8,cipher_A15B12;
    cipher_matrix_jiang cipher_A12B1,cipher_A13B5,cipher_A14B9,cipher_A15B13;
    cipher_matrix_jiang cipher_A12B2,cipher_A13B6,cipher_A14B10,cipher_A15B14;
    cipher_matrix_jiang cipher_A12B3,cipher_A13B7,cipher_A14B11,cipher_A15B15;

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A00,cipher_B00,cipher_A0B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A00,cipher_B1_1,cipher_A0B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A00,cipher_B22,cipher_A0B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A00,cipher_B33,cipher_A0B3);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A44,cipher_B00,cipher_A4B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A44,cipher_B1_1,cipher_A4B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A44,cipher_B22,cipher_A4B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A44,cipher_B33,cipher_A4B3);
        
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A88,cipher_B00,cipher_A8B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A88,cipher_B1_1,cipher_A8B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A88,cipher_B22,cipher_A8B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A88,cipher_B33,cipher_A8B3);
        
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1212,cipher_B00,cipher_A12B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1212,cipher_B1_1,cipher_A12B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1212,cipher_B22,cipher_A12B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1212,cipher_B33,cipher_A12B3);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1_1,cipher_B44,cipher_A1B4);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1_1,cipher_B55,cipher_A1B5);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1_1,cipher_B66,cipher_A1B6);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1_1,cipher_B77,cipher_A1B7);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A55,cipher_B44,cipher_A5B4);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A55,cipher_B55,cipher_A5B5);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A55,cipher_B66,cipher_A5B6);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A55,cipher_B77,cipher_A5B7);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A99,cipher_B44,cipher_A9B4);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A99,cipher_B55,cipher_A9B5);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A99,cipher_B66,cipher_A9B6);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A99,cipher_B77,cipher_A9B7);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1313,cipher_B44,cipher_A13B4);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1313,cipher_B55,cipher_A13B5);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1313,cipher_B66,cipher_A13B6);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1313,cipher_B77,cipher_A13B7);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A22,cipher_B88,cipher_A2B8);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A22,cipher_B99,cipher_A2B9);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A22,cipher_B1010,cipher_A2B10);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A22,cipher_B1111,cipher_A2B11);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A66,cipher_B88,cipher_A6B8);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A66,cipher_B99,cipher_A6B9);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A66,cipher_B1010,cipher_A6B10);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A66,cipher_B1111,cipher_A6B11);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1010,cipher_B88,cipher_A10B8);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1010,cipher_B99,cipher_A10B9);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1010,cipher_B1010,cipher_A10B10);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1010,cipher_B1111,cipher_A10B11);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1414,cipher_B88,cipher_A14B8);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1414,cipher_B99,cipher_A14B9);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1414,cipher_B1010,cipher_A14B10);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1414,cipher_B1111,cipher_A14B11);
       
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A33,cipher_B1212,cipher_A3B12);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A33,cipher_B1313,cipher_A3B13);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A33,cipher_B1414,cipher_A3B14);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A33,cipher_B1515,cipher_A3B15);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A77,cipher_B1212,cipher_A7B12);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A77,cipher_B1313,cipher_A7B13);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A77,cipher_B1414,cipher_A7B14);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A77,cipher_B1515,cipher_A7B15);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1111,cipher_B1212,cipher_A11B12);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1111,cipher_B1313,cipher_A11B13);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1111,cipher_B1414,cipher_A11B14);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1111,cipher_B1515,cipher_A11B15);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1515,cipher_B1212,cipher_A15B12);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1515,cipher_B1313,cipher_A15B13);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1515,cipher_B1414,cipher_A15B14);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1515,cipher_B1515,cipher_A15B15);
   
    seal::Ciphertext cipher_A0B0_ciphertext,cipher_A1B4_ciphertext,cipher_A2B8_ciphertext,cipher_A3B12_ciphertext;
    seal::Ciphertext cipher_A0B1_ciphertext,cipher_A1B5_ciphertext,cipher_A2B9_ciphertext,cipher_A3B13_ciphertext;
    seal::Ciphertext cipher_A0B2_ciphertext,cipher_A1B6_ciphertext,cipher_A2B10_ciphertext,cipher_A3B14_ciphertext;
    seal::Ciphertext cipher_A0B3_ciphertext,cipher_A1B7_ciphertext,cipher_A2B11_ciphertext,cipher_A3B15_ciphertext;
    seal::Ciphertext cipher_A4B0_ciphertext,cipher_A5B4_ciphertext,cipher_A6B8_ciphertext,cipher_A7B12_ciphertext;
    seal::Ciphertext cipher_A4B1_ciphertext,cipher_A5B5_ciphertext,cipher_A6B9_ciphertext,cipher_A7B13_ciphertext;
    seal::Ciphertext cipher_A4B2_ciphertext,cipher_A5B6_ciphertext,cipher_A6B10_ciphertext,cipher_A7B14_ciphertext;
    seal::Ciphertext cipher_A4B3_ciphertext,cipher_A5B7_ciphertext,cipher_A6B11_ciphertext,cipher_A7B15_ciphertext;

    seal::Ciphertext cipher_A8B0_ciphertext,cipher_A9B4_ciphertext,cipher_A10B8_ciphertext,cipher_A11B12_ciphertext;
    seal::Ciphertext cipher_A8B1_ciphertext,cipher_A9B5_ciphertext,cipher_A10B9_ciphertext,cipher_A11B13_ciphertext;
    seal::Ciphertext cipher_A8B2_ciphertext,cipher_A9B6_ciphertext,cipher_A10B10_ciphertext,cipher_A11B14_ciphertext;
    seal::Ciphertext cipher_A8B3_ciphertext,cipher_A9B7_ciphertext,cipher_A10B11_ciphertext,cipher_A11B15_ciphertext;

    seal::Ciphertext cipher_A12B0_ciphertext,cipher_A13B4_ciphertext,cipher_A14B8_ciphertext,cipher_A15B12_ciphertext;
    seal::Ciphertext cipher_A12B1_ciphertext,cipher_A13B5_ciphertext,cipher_A14B9_ciphertext,cipher_A15B13_ciphertext;
    seal::Ciphertext cipher_A12B2_ciphertext,cipher_A13B6_ciphertext,cipher_A14B10_ciphertext,cipher_A15B14_ciphertext;
    seal::Ciphertext cipher_A12B3_ciphertext,cipher_A13B7_ciphertext,cipher_A14B11_ciphertext,cipher_A15B15_ciphertext;
    
    seal::Ciphertext cipher_temp_0,cipher_temp_1,cipher_temp_2,cipher_temp_3;
    seal::Ciphertext cipher_temp_4,cipher_temp_5,cipher_temp_6,cipher_temp_7;
    seal::Ciphertext cipher_temp_8,cipher_temp_9,cipher_temp_10,cipher_temp_11;
    seal::Ciphertext cipher_temp_12,cipher_temp_13,cipher_temp_14,cipher_temp_15;

    cipher_A0B0_ciphertext=cipher_A0B0.get_cipher_matrix();
    cipher_A0B1_ciphertext=cipher_A0B1.get_cipher_matrix();
    cipher_A0B2_ciphertext=cipher_A0B2.get_cipher_matrix();
    cipher_A0B3_ciphertext=cipher_A0B3.get_cipher_matrix();
    cipher_A4B0_ciphertext=cipher_A4B0.get_cipher_matrix();
    cipher_A4B1_ciphertext=cipher_A4B1.get_cipher_matrix();
    cipher_A4B2_ciphertext=cipher_A4B2.get_cipher_matrix();
    cipher_A4B3_ciphertext=cipher_A4B3.get_cipher_matrix();
    cipher_A8B0_ciphertext=cipher_A8B0.get_cipher_matrix();
    cipher_A8B1_ciphertext=cipher_A8B1.get_cipher_matrix();
    cipher_A8B2_ciphertext=cipher_A8B2.get_cipher_matrix();
    cipher_A8B3_ciphertext=cipher_A8B3.get_cipher_matrix();
    cipher_A12B0_ciphertext=cipher_A12B0.get_cipher_matrix();
    cipher_A12B1_ciphertext=cipher_A12B1.get_cipher_matrix();
    cipher_A12B2_ciphertext=cipher_A12B2.get_cipher_matrix();
    cipher_A12B3_ciphertext=cipher_A12B3.get_cipher_matrix();

    cipher_A1B4_ciphertext=cipher_A1B4.get_cipher_matrix();
    cipher_A1B5_ciphertext=cipher_A1B5.get_cipher_matrix();
    cipher_A1B6_ciphertext=cipher_A1B6.get_cipher_matrix();
    cipher_A1B7_ciphertext=cipher_A1B7.get_cipher_matrix();
    cipher_A5B4_ciphertext=cipher_A5B4.get_cipher_matrix();
    cipher_A5B5_ciphertext=cipher_A5B5.get_cipher_matrix();
    cipher_A5B6_ciphertext=cipher_A5B6.get_cipher_matrix();
    cipher_A5B7_ciphertext=cipher_A5B7.get_cipher_matrix();
    cipher_A9B4_ciphertext=cipher_A9B4.get_cipher_matrix();
    cipher_A9B5_ciphertext=cipher_A9B5.get_cipher_matrix();
    cipher_A9B6_ciphertext=cipher_A9B6.get_cipher_matrix();
    cipher_A9B7_ciphertext=cipher_A9B7.get_cipher_matrix();
    cipher_A13B4_ciphertext=cipher_A13B4.get_cipher_matrix();
    cipher_A13B5_ciphertext=cipher_A13B5.get_cipher_matrix();
    cipher_A13B6_ciphertext=cipher_A13B6.get_cipher_matrix();
    cipher_A13B7_ciphertext=cipher_A13B7.get_cipher_matrix();

    cipher_A2B8_ciphertext=cipher_A2B8.get_cipher_matrix();
    cipher_A2B9_ciphertext=cipher_A2B9.get_cipher_matrix();
    cipher_A2B10_ciphertext=cipher_A2B10.get_cipher_matrix();
    cipher_A2B11_ciphertext=cipher_A2B11.get_cipher_matrix();
    cipher_A6B8_ciphertext=cipher_A6B8.get_cipher_matrix();
    cipher_A6B9_ciphertext=cipher_A6B9.get_cipher_matrix();
    cipher_A6B10_ciphertext=cipher_A6B10.get_cipher_matrix();
    cipher_A6B11_ciphertext=cipher_A6B11.get_cipher_matrix();
    cipher_A10B8_ciphertext=cipher_A10B8.get_cipher_matrix();
    cipher_A10B9_ciphertext=cipher_A10B9.get_cipher_matrix();
    cipher_A10B10_ciphertext=cipher_A10B10.get_cipher_matrix();
    cipher_A10B11_ciphertext=cipher_A10B11.get_cipher_matrix();
    cipher_A14B8_ciphertext=cipher_A14B8.get_cipher_matrix();
    cipher_A14B9_ciphertext=cipher_A14B9.get_cipher_matrix();
    cipher_A14B10_ciphertext=cipher_A14B10.get_cipher_matrix();
    cipher_A14B11_ciphertext=cipher_A14B11.get_cipher_matrix();

    cipher_A3B12_ciphertext=cipher_A3B12.get_cipher_matrix();
    cipher_A3B13_ciphertext=cipher_A3B13.get_cipher_matrix();
    cipher_A3B14_ciphertext=cipher_A3B14.get_cipher_matrix();
    cipher_A3B15_ciphertext=cipher_A3B15.get_cipher_matrix();
    cipher_A7B12_ciphertext=cipher_A7B12.get_cipher_matrix();
    cipher_A7B13_ciphertext=cipher_A7B13.get_cipher_matrix();
    cipher_A7B14_ciphertext=cipher_A7B14.get_cipher_matrix();
    cipher_A7B15_ciphertext=cipher_A7B15.get_cipher_matrix();
    cipher_A11B12_ciphertext=cipher_A11B12.get_cipher_matrix();
    cipher_A11B13_ciphertext=cipher_A11B13.get_cipher_matrix();
    cipher_A11B14_ciphertext=cipher_A11B14.get_cipher_matrix();
    cipher_A11B15_ciphertext=cipher_A11B15.get_cipher_matrix();
    cipher_A15B12_ciphertext=cipher_A15B12.get_cipher_matrix();
    cipher_A15B13_ciphertext=cipher_A15B13.get_cipher_matrix();
    cipher_A15B14_ciphertext=cipher_A15B14.get_cipher_matrix();
    cipher_A15B15_ciphertext=cipher_A15B15.get_cipher_matrix();

    evaluator.add(cipher_A0B0_ciphertext,cipher_A1B4_ciphertext,cipher_temp_0);
    evaluator.add_inplace(cipher_temp_0,cipher_A2B8_ciphertext);
    evaluator.add_inplace(cipher_temp_0,cipher_A3B12_ciphertext);
    evaluator.add(cipher_A0B1_ciphertext,cipher_A1B5_ciphertext,cipher_temp_1);
    evaluator.add_inplace(cipher_temp_1,cipher_A2B9_ciphertext);
    evaluator.add_inplace(cipher_temp_1,cipher_A3B13_ciphertext);
    evaluator.add(cipher_A0B2_ciphertext,cipher_A1B6_ciphertext,cipher_temp_2);
    evaluator.add_inplace(cipher_temp_2,cipher_A2B10_ciphertext);
    evaluator.add_inplace(cipher_temp_2,cipher_A3B14_ciphertext);
    evaluator.add(cipher_A0B3_ciphertext,cipher_A1B7_ciphertext,cipher_temp_3);
    evaluator.add_inplace(cipher_temp_3,cipher_A2B11_ciphertext);
    evaluator.add_inplace(cipher_temp_3,cipher_A3B15_ciphertext);

    evaluator.add(cipher_A4B0_ciphertext,cipher_A5B4_ciphertext,cipher_temp_4);
    evaluator.add_inplace(cipher_temp_4,cipher_A6B8_ciphertext);
    evaluator.add_inplace(cipher_temp_4,cipher_A7B12_ciphertext);
    evaluator.add(cipher_A4B1_ciphertext,cipher_A5B5_ciphertext,cipher_temp_5);
    evaluator.add_inplace(cipher_temp_5,cipher_A6B9_ciphertext);
    evaluator.add_inplace(cipher_temp_5,cipher_A7B13_ciphertext);
    evaluator.add(cipher_A4B2_ciphertext,cipher_A5B6_ciphertext,cipher_temp_6);
    evaluator.add_inplace(cipher_temp_6,cipher_A6B10_ciphertext);
    evaluator.add_inplace(cipher_temp_6,cipher_A7B14_ciphertext);
    evaluator.add(cipher_A4B3_ciphertext,cipher_A5B7_ciphertext,cipher_temp_7);
    evaluator.add_inplace(cipher_temp_7,cipher_A6B11_ciphertext);
    evaluator.add_inplace(cipher_temp_7,cipher_A7B15_ciphertext);

    evaluator.add(cipher_A8B0_ciphertext,cipher_A9B4_ciphertext,cipher_temp_8);
    evaluator.add_inplace(cipher_temp_8,cipher_A10B8_ciphertext);
    evaluator.add_inplace(cipher_temp_8,cipher_A11B12_ciphertext);
    evaluator.add(cipher_A8B1_ciphertext,cipher_A9B5_ciphertext,cipher_temp_9);
    evaluator.add_inplace(cipher_temp_9,cipher_A10B9_ciphertext);
    evaluator.add_inplace(cipher_temp_9,cipher_A11B13_ciphertext);
    evaluator.add(cipher_A8B2_ciphertext,cipher_A9B6_ciphertext,cipher_temp_10);
    evaluator.add_inplace(cipher_temp_10,cipher_A10B10_ciphertext);
    evaluator.add_inplace(cipher_temp_10,cipher_A11B14_ciphertext);
    evaluator.add(cipher_A8B3_ciphertext,cipher_A9B7_ciphertext,cipher_temp_11);
    evaluator.add_inplace(cipher_temp_11,cipher_A10B11_ciphertext);
    evaluator.add_inplace(cipher_temp_11,cipher_A11B15_ciphertext);

    evaluator.add(cipher_A12B0_ciphertext,cipher_A13B4_ciphertext,cipher_temp_12);
    evaluator.add_inplace(cipher_temp_12,cipher_A14B8_ciphertext);
    evaluator.add_inplace(cipher_temp_12,cipher_A15B12_ciphertext);
    evaluator.add(cipher_A12B1_ciphertext,cipher_A13B5_ciphertext,cipher_temp_13);
    evaluator.add_inplace(cipher_temp_13,cipher_A14B9_ciphertext);
    evaluator.add_inplace(cipher_temp_13,cipher_A15B13_ciphertext);
    evaluator.add(cipher_A12B2_ciphertext,cipher_A13B6_ciphertext,cipher_temp_14);
    evaluator.add_inplace(cipher_temp_14,cipher_A14B10_ciphertext);
    evaluator.add_inplace(cipher_temp_14,cipher_A15B14_ciphertext);
    evaluator.add(cipher_A12B3_ciphertext,cipher_A13B7_ciphertext,cipher_temp_15);
    evaluator.add_inplace(cipher_temp_15,cipher_A14B11_ciphertext);
    evaluator.add_inplace(cipher_temp_15,cipher_A15B15_ciphertext);
    
    
    //compute A^2+A
    seal::Ciphertext cipher_temp_0_A0,cipher_temp_1_A1,cipher_temp_2_A2,cipher_temp_3_A3;
    seal::Ciphertext cipher_temp_4_A4,cipher_temp_5_A5,cipher_temp_6_A6,cipher_temp_7_A7;
    seal::Ciphertext cipher_temp_8_A8,cipher_temp_9_A9,cipher_temp_10_A10,cipher_temp_11_A11;
    seal::Ciphertext cipher_temp_12_A12,cipher_temp_13_A13,cipher_temp_14_A14,cipher_temp_15_A15;

    match_levels(context,cipher_A0_ciphertext,cipher_temp_0,evaluator);
    evaluator.add(cipher_temp_0,cipher_A0_ciphertext,cipher_temp_0_A0);
    match_levels(context,cipher_A1_ciphertext,cipher_temp_1,evaluator);
    evaluator.add(cipher_temp_1,cipher_A1_ciphertext,cipher_temp_1_A1);
    match_levels(context,cipher_A2_ciphertext,cipher_temp_2,evaluator);
    evaluator.add(cipher_temp_2,cipher_A2_ciphertext,cipher_temp_2_A2);
    match_levels(context,cipher_A3_ciphertext,cipher_temp_3,evaluator);
    evaluator.add(cipher_temp_3,cipher_A3_ciphertext,cipher_temp_3_A3);

    match_levels(context,cipher_A4_ciphertext,cipher_temp_4,evaluator);
    evaluator.add(cipher_temp_4,cipher_A4_ciphertext,cipher_temp_4_A4);
    match_levels(context,cipher_A5_ciphertext,cipher_temp_5,evaluator);
    evaluator.add(cipher_temp_5,cipher_A5_ciphertext,cipher_temp_5_A5);
    match_levels(context,cipher_A6_ciphertext,cipher_temp_6,evaluator);
    evaluator.add(cipher_temp_6,cipher_A6_ciphertext,cipher_temp_6_A6);
    match_levels(context,cipher_A7_ciphertext,cipher_temp_7,evaluator);
    evaluator.add(cipher_temp_7,cipher_A7_ciphertext,cipher_temp_7_A7);

    match_levels(context,cipher_A8_ciphertext,cipher_temp_8,evaluator);
    evaluator.add(cipher_temp_8,cipher_A8_ciphertext,cipher_temp_8_A8);
    match_levels(context,cipher_A9_ciphertext,cipher_temp_9,evaluator);
    evaluator.add(cipher_temp_9,cipher_A9_ciphertext,cipher_temp_9_A9);
    match_levels(context,cipher_A10_ciphertext,cipher_temp_10,evaluator);
    evaluator.add(cipher_temp_10,cipher_A10_ciphertext,cipher_temp_10_A10);
    match_levels(context,cipher_A11_ciphertext,cipher_temp_11,evaluator);
    evaluator.add(cipher_temp_11,cipher_A11_ciphertext,cipher_temp_11_A11);

    match_levels(context,cipher_A12_ciphertext,cipher_temp_4,evaluator);
    evaluator.add(cipher_temp_12,cipher_A12_ciphertext,cipher_temp_12_A12);
    match_levels(context,cipher_A13_ciphertext,cipher_temp_13,evaluator);
    evaluator.add(cipher_temp_13,cipher_A13_ciphertext,cipher_temp_13_A13);
    match_levels(context,cipher_A14_ciphertext,cipher_temp_14,evaluator);
    evaluator.add(cipher_temp_14,cipher_A14_ciphertext,cipher_temp_14_A14);
    match_levels(context,cipher_A15_ciphertext,cipher_temp_15,evaluator);
    evaluator.add(cipher_temp_15,cipher_A15_ciphertext,cipher_temp_15_A15);
    destination.push_back(cipher_temp_0_A0);
    destination.push_back(cipher_temp_1_A1);
    destination.push_back(cipher_temp_2_A2);
    destination.push_back(cipher_temp_3_A3);
    destination.push_back(cipher_temp_4_A4);
    destination.push_back(cipher_temp_5_A5);
    destination.push_back(cipher_temp_6_A6);
    destination.push_back(cipher_temp_7_A7);
    destination.push_back(cipher_temp_8_A8);
    destination.push_back(cipher_temp_9_A9);
    destination.push_back(cipher_temp_10_A10);
    destination.push_back(cipher_temp_11_A11);
    destination.push_back(cipher_temp_12_A12);
    destination.push_back(cipher_temp_13_A13);
    destination.push_back(cipher_temp_14_A14);
    destination.push_back(cipher_temp_15_A15);
    return destination;
}                          