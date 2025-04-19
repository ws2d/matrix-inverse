
#include "matrix.h"
#include <seal/seal.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <seal/util/numth.h>


// using namespace seal;
using namespace seal::util;
using namespace std;

void eye_matrix_generator(std::size_t n, std::size_t m, Matrix<double> & A){
    A.resize(n,m);
    for (size_t i = 0; i < m; i++){
        for(size_t j = 0; j < n; j++){
            A.set(j, i, 0);
            if(i==j){
                A.set(j,i,1);
            }
        }
    }
}

void random_matrix_generator(std::size_t n, std::size_t m, Matrix<double> & A){
    A.resize(n, m);
    for (size_t i = 0; i < m; i++){
        for(size_t j = 0; j < n; j++){
            A.set(j, i, pow(-1, i+j)*rand()/(10*(pow(2, 35))));
            if(i==j){
                A.set(j,i,10);
            }
        }
    }
    A.print();
    cout<<"wssssss"<<endl;

    // // 将矩阵存储到文件
    // std::ofstream outFile("matrix_128.txt");
    // for (size_t i = 0; i < m; i++) {
    //     for (size_t j = 0; j < n; j++) {
    //         outFile << A.get(j, i) << (j < n - 1 ? " " : "");
    //     }
    //     outFile << std::endl;
    // }

    // outFile.close();
    // cout<<"打印矩阵A:"<<endl;
    // A.print();
}

void random_tae_matrix_generator(std::size_t n, std::size_t m, Matrix<double> & A){
    A.resize(n, m);
    for (size_t i = 0; i < m; i++){
        for(size_t j = 0; j < n; j++){
            A.set(j, i, (pow(-1, i + j) * (0.01 + (0.1 - 0.01) * static_cast<double>(rand())/20) / RAND_MAX));
            if(i==j){
                A.set(j,i,0.9);
            }
        }
    }

}
void random_block_matrix_generator(std::size_t n, std::size_t m, Matrix<double> & A){
    A.resize(n, m);
    for (size_t i = 0; i < m; i++){
        for(size_t j = 0; j < n; j++){
            A.set(j, i, (pow(-1, i + j) * (0.01 + (0.1 - 0.01) * static_cast<double>(rand())/400) / RAND_MAX));
            if(i==j){
                A.set(j,i,0.00001);
            }
        }
    }
    A.print();
    cout<<"wssssss"<<endl;

    // 将矩阵存储到文件
    std::ofstream outFile("matrix_256.txt");
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            outFile << A.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile << std::endl;
    }

    outFile.close();
    cout<<"打印矩阵A:"<<endl;
    A.print();
}

void random_block_matrix_generator_512(std::size_t n, std::size_t m, Matrix<double> & A){
    A.resize(n, m);
    for (size_t i = 0; i < m; i++){
        for(size_t j = 0; j < n; j++){
            A.set(j, i, (pow(-1, i + j) * (0.01 + (0.1 - 0.01) * static_cast<double>(rand())/300) / RAND_MAX));
            if(i==j){
                A.set(j,i,0.00001);
            }
        }
    }
    A.print();
    cout<<"wssssss"<<endl;

    // 将矩阵存储到文件
    std::ofstream outFile("matrix_512.txt");
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            outFile << A.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile << std::endl;
    }

    outFile.close();
    cout<<"打印矩阵A:"<<endl;
    A.print();
}
void random_pack_matrix_generator_new_encoding(std::size_t n, std::size_t m, Matrix<double> & A){
    A.resize(n, m);
    // size_t p=1;

    for (size_t i = 0; i < m; i++){
        for(size_t j = 0; j < n; j++){
            // int random_value = rand() % 10 + 1;
            // A.set(i, j, random_value);
            // A.set(i, j, p);
            // p++;
            // A.set(j, i, pow(-1, i+j)*rand()/(pow(2, 35)));
            A.set(j, i, (pow(-1, i + j) * (0.01 + (0.1 - 0.01) * static_cast<double>(rand())/10) / RAND_MAX));
            // d=4,8,16, 2,0.99
            //d=32 64, 5,0.99
            // d=,128, 为10,0.999
            //加减10^-n，n取为数量级
            if(i==j){
                A.set(j,i,0.999);
            }
        }
    }
}
void random_pack_triangle_matrix_generator(std::size_t n, std::size_t m, Matrix<double> & A){
    A.resize(n, m);


    for (size_t i = 0; i < m; i++){
        for(size_t j = 0; j < n; j++){
            // int random_value = rand() % 10 + 1;
            // A.set(i, j, random_value);
            if(j<i){
                // int random_value = rand() % 10 + 1;
                double random_value = 1.0 + static_cast<double>(rand()) /(RAND_MAX /(10.0-1.0));
                // 创建一个 stringstream 对象
                std::stringstream ss;
                // 设置小数点后保留3位
                ss << std::fixed << std::setprecision(3) << random_value;

                // 将格式化后的字符串转换回 double
                double formatted_value = std::stod(ss.str());
                A.set(i, j, formatted_value);
            }


        }
    }
}

void random_pack_matrix_generator(std::size_t n, std::size_t m, Matrix<double> & A){
    A.resize(n, m);
    // size_t p=1;

    for (size_t i = 0; i < m; i++){
        for(size_t j = 0; j < n; j++){
            // int random_value = rand() % 10 + 1;
            // A.set(i, j, random_value);
            // A.set(i, j, p);
            // p++;
            A.set(j, i, pow(-1, i+j)*rand()/(pow(2, 35)));
            if(i==j){
                A.set(j,i,0.8);
            }
        }
    }
    // A.print();
    // cout<<"pack_example"<<endl;
}

void store_matrix(std::size_t n, std::size_t m, Matrix<double> & A){
    std::ofstream outFile("matrix_result_512.txt");
     // 设置输出格式为小数形式，并设置精度
    outFile << std::fixed << std::setprecision(12); // 可以根据需要调整精度
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            outFile << A.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile << std::endl;
    }

    outFile.close();
    cout<<"打印矩阵A:"<<endl;
}

void store_vector(std::vector<double>& vec)
{
    // 定义文件路径和文件名
    std::string filePath = "../build/dia_triangle/result_triangle_(4,4).txt";  // 修改为你的目标路径

    // 打开文件
    std::ofstream outFile(filePath);

    // 检查文件是否成功打开
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filePath << std::endl;
        return;
    }

    // 将向量内容写入文件，用空格分隔每个值
    for (size_t i = 0; i < vec.size(); ++i) {
        outFile << vec[i];
        if (i < vec.size() - 1) {
            outFile << " ";  // 在每个值之间添加空格
        }
    }

    // 关闭文件
    outFile.close();
    // std::cout << "Vector saved to file: " << filePath << std::endl;
}
void store_vector_triangle(std::vector<double>& vec)
{
    // 定义文件路径和文件名
    std::string filePath = "./triangle_input_16384.txt";  // 修改为你的目标路径

    // 打开文件
    std::ofstream outFile(filePath);

    // 检查文件是否成功打开
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filePath << std::endl;
        return;
    }

    // 将向量内容写入文件，用空格分隔每个值
    for (size_t i = 0; i < vec.size(); ++i) {
        outFile << vec[i];
        if (i < vec.size() - 1) {
            outFile << " ";  // 在每个值之间添加空格
        }
    }

    // 关闭文件
    outFile.close();
    // std::cout << "Vector saved to file: " << filePath << std::endl;
}

void store_vector_new_encoding(std::vector<double>& vec)
{
    // 定义文件路径和文件名
    std::string filePath = "./tae_exp/encoding_input(128,128)";  // 修改为你的目标路径

    // 打开文件
    std::ofstream outFile(filePath);

    // 检查文件是否成功打开
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filePath << std::endl;
        return;
    }

    // 将向量内容写入文件，用空格分隔每个值
    for (size_t i = 0; i < vec.size(); ++i) {
        outFile << vec[i];
        if (i < vec.size() - 1) {
            outFile << " ";  // 在每个值之间添加空格
        }
    }

    // 关闭文件
    outFile.close();
    // std::cout << "Vector saved to file: " << filePath << std::endl;
}

//读取文件
std::vector<double> read_data_from_file(const std::string& file_path) {
    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << file_path << std::endl;
        throw std::runtime_error("无法打开文件");
    }

    std::vector<double> data;
    double value;

    // 逐个读取文件中的数值
    while (file >> value) {
        data.push_back(value);
    }

    file.close();
    return data;
}

void store_vector_new_encoding_over(std::vector<double>& vec)
{
    // 定义文件路径和文件名
    std::string filePath = "../build/dia_triangle/result_tae_(4,4).txt";  // 修改为你的目标路径

    // 打开文件
    std::ofstream outFile(filePath);

    // 检查文件是否成功打开
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to open file " << filePath << std::endl;
        return;
    }
    // 设置输出格式：固定小数点格式，并保留足够的小数位数
    outFile << std::fixed << std::setprecision(10); // 保留 10 位小数
    // 将向量内容写入文件，用空格分隔每个值
    for (size_t i = 0; i < vec.size(); ++i) {
        outFile << vec[i];
        if (i < vec.size() - 1) {
            outFile << " ";  // 在每个值之间添加空格
        }
    }

    // 关闭文件
    outFile.close();
    // std::cout << "Vector saved to file: " << filePath << std::endl;
}


//pack

// 解密函数
void decrypt_and_print(seal::Ciphertext& cipher, 
                                      seal::Decryptor& decryptor, 
                                      seal::CKKSEncoder& encoder) {
    // 创建明文对象用于存储解密结果
    seal::Plaintext plain_over;

    // 解密密文
    decryptor.decrypt(cipher, plain_over);

    // 解码明文为向量
    std::vector<double> result;
    encoder.decode(plain_over, result);
    // store_vector(result);

    // 打印解密后的向量前16位
    std::cout << "解密后: 向量前16位" << std::endl;
    for (size_t i = 0; i < result.size() && i < 32; ++i) {
        if(i%16==0&&i!=0)
        {
            cout<<endl;
        }
        std::cout <<std::fixed<<std::setprecision(10)<< result[i] << " ";
    }
    cout<<endl;
}

// 解密函数
void decrypt_and_print_new_encoding(seal::Ciphertext& cipher, 
                                      seal::Decryptor& decryptor, 
                                      seal::CKKSEncoder& encoder,
                                      size_t d) {
    size_t step=encoder.slot_count()/d/d;
    // 创建明文对象用于存储解密结果
    seal::Plaintext plain_over;

    // 解密密文
    decryptor.decrypt(cipher, plain_over);

    // 解码明文为向量
    std::vector<double> result;
    encoder.decode(plain_over, result);

    // store_vector_new_encoding_over(result);

    // 打印解密后的向量前16位
    std::cout << "解密后: 每块第1位:" << std::endl;
    for (size_t i = 0; i < result.size(); ++i) {
        if(i%step==0)
        {
            std::cout << result[i] << " ";
        }
    }
    cout<<endl;
}
void decrypt_and_print_new_encoding_over(seal::Ciphertext& cipher, 
                                      seal::Decryptor& decryptor, 
                                      seal::CKKSEncoder& encoder,
                                      std::vector<double> result_I,
                                      size_t d) {
    size_t step=encoder.slot_count()/d/d;
    // 创建明文对象用于存储解密结果
    seal::Plaintext plain_over;

    // 解密密文
    decryptor.decrypt(cipher, plain_over);

    // 解码明文为向量
    std::vector<double> result;
    encoder.decode(plain_over, result);
    for (size_t i = 0; i < result_I.size(); ++i) {
        result[i]=result_I[i]+result[i];
    }
    // store_vector_new_encoding_over(result);

    // 打印解密后的向量前16位
    std::cout << "解密后: 每块第1位:" << std::endl;
    for (size_t i = 0; i < result.size(); ++i) {
        if(i%step==0)
        {
            std::cout << result[i] << " ";
        }
    }
    cout<<endl;
}


seal::Ciphertext rotate_as_16(seal::Ciphertext& cipher_matrix,int d, double scale, 
    seal::CKKSEncoder &encoder, 
    seal::Evaluator &evaluator, 
    seal::Decryptor &decryptor,
    seal::GaloisKeys &galois_keys, 
    seal::RelinKeys &relin_keys)
{
    seal::Ciphertext cipher_tmp,cipher_tmp1,cipher_tmp_left,cipher_tmp_right,rotate_result;
    seal::Plaintext plain_tmp,plain_tmp1;
    vector<double> vec_tmp(16,0);
    size_t step=encoder.slot_count()/16;
    // cout<<"step"<<step<<endl;
    if(d==0){
        vector<double> rotate_all(16384, 1);
        seal::Plaintext plain_tmp;
        encoder.encode(rotate_all, cipher_matrix.scale(), plain_tmp);
        evaluator.mod_switch_to_inplace(plain_tmp, cipher_matrix.parms_id());
        evaluator.multiply_plain(cipher_matrix, plain_tmp, rotate_result);
        evaluator.relinearize_inplace(rotate_result, relin_keys);
        evaluator.rescale_to_next_inplace(rotate_result);



        // rotate_result=cipher_matrix;
        // evaluator.mod_switch_to_next_inplace(rotate_result);

        return rotate_result;
    }
    else{

        if(d<0){
            d=d+16;
        }
        for(int i=0 ;i<d;i++){
            vec_tmp[i]=1;
        }
        vector<double> result;
        // 重复 1024 次
        for (size_t s = 0; s < step; s++) {
            result.insert(result.end(), vec_tmp.begin(), vec_tmp.end());
        }
        // 将结果保存回 vec_tmp
        vec_tmp = result;
        // for (int s = 0; s < 32; s++) {
        //     cout<<vec_tmp[s]<<" ";
        // }

        encoder.encode(vec_tmp, cipher_matrix.scale(), plain_tmp);
        evaluator.mod_switch_to_inplace(plain_tmp, cipher_matrix.parms_id());
        evaluator.multiply_plain(cipher_matrix, plain_tmp, cipher_tmp);
        evaluator.rotate_vector(cipher_tmp, d-16, galois_keys, cipher_tmp_left);
        // decrypt_and_print(cipher_tmp_left, decryptor, encoder);
        //右旋开始
        vector<double> vec_tmp1(16,0);
        for(int i=0;i<16-d;i++)
        {
            vec_tmp1[i]=1;
        }
        vector<double> result1;
        // 重复 1024 次
        for (size_t s = 0; s < step; s++) {
            result1.insert(result1.end(), vec_tmp1.begin(), vec_tmp1.end());
        }
        // 将结果保存回 vec_tmp
        vec_tmp1 = result1;

        evaluator.rotate_vector(cipher_matrix, d, galois_keys, cipher_tmp1);
        encoder.encode(vec_tmp1, cipher_matrix.scale(), plain_tmp1);
        evaluator.mod_switch_to_inplace(plain_tmp1, cipher_matrix.parms_id());
        evaluator.multiply_plain(cipher_tmp1, plain_tmp1, cipher_tmp_right);
        evaluator.add(cipher_tmp_left,cipher_tmp_right,rotate_result);
        evaluator.relinearize_inplace(rotate_result, relin_keys);
        evaluator.rescale_to_next_inplace(rotate_result);
        // cipher_matrix=rotate_result;


        return rotate_result;

    }
    return cipher_matrix;
}

seal::Ciphertext pack_change_cipher_to_matrix_A(
    seal::SEALContext context,
    seal::Ciphertext& cipher_matrix, 
    double scale, 
    seal::CKKSEncoder &encoder, 
    seal::Evaluator &evaluator, 
    seal::Decryptor &decryptor,
    seal::GaloisKeys &galois_keys, 
    seal::RelinKeys &relin_keys)
{
    size_t d=4;
    Matrix<double> u_sigma;
    u_sigma.generate_u_sigma(d, d);
    seal::Ciphertext cipher_tmp,cipher_tmp_1;//save intermediate variables.
    vector<seal::Ciphertext> cipher_rotate,cipher_result;//save rotate result vector
    seal::Plaintext plain_tmp;
    vector<double> vec_tmp;
    seal::Ciphertext destination;

    //step1. determining the dimensions
    int rotate_size = 2 * d - 1;
    int rotate_outside_size = ceil(sqrt(rotate_size));
    int rotate_inside_size = ceil(double(rotate_size) / double(rotate_outside_size));
    // cout<<rotate_inside_size<<rotate_outside_size<<endl;

    // step2. pre-rotate ciphertext
    

    for (int i = 0; i < rotate_inside_size; i++) {

        cipher_tmp= rotate_as_16(cipher_matrix,-d+d*d+i+1, scale, encoder, evaluator, decryptor,galois_keys, relin_keys);
        // cout<<-d+d*d+i+1<<endl;
        // decrypt_and_print(cipher_tmp, decryptor, encoder);
        cipher_rotate.push_back(cipher_tmp);
    }

    //step3. starting step function 
    for (int i = 0; i < rotate_outside_size; i++) {
        vector<seal::Ciphertext> cipher_vector;
        for (int j = 0; j < rotate_inside_size; j++) {
            vec_tmp = u_sigma.diag_vector(i * rotate_inside_size + j-d+1,d*d);
            // for(size_t s=0;s<vec_tmp.size();s++)
            // {
            //     cout<<vec_tmp[s]<<" ";
            // }
            // cout<<endl;
            // cout<<"vec_tmp.size:"<<vec_tmp.size()<<endl;
            if (std::all_of(vec_tmp.begin(), vec_tmp.end(), [](double num) { return num == 0; }) == 1) {
                continue;
            }
            else {
                std::rotate(vec_tmp.rbegin(), vec_tmp.rbegin() + i * rotate_inside_size, vec_tmp.rend());
                // for(size_t s=0;s<vec_tmp.size();s++)
                // {
                //     cout<<vec_tmp[s]<<" ";
                // }
                // cout<<endl;
                vector<double> result;
            
                // 重复 1024 次
                for (int s = 0; s < 1024; s++) {
                    result.insert(result.end(), vec_tmp.begin(), vec_tmp.end());
                }
                // 将结果保存回 vec_tmp
                vec_tmp = result;
                encoder.encode(vec_tmp, cipher_rotate[0].scale(), plain_tmp);
                evaluator.mod_switch_to_inplace(plain_tmp, cipher_rotate[0].parms_id());
                evaluator.multiply_plain(cipher_rotate[j], plain_tmp, cipher_tmp);
                
                cipher_vector.push_back(cipher_tmp);
                
            }
        }
        evaluator.add_many(cipher_vector, cipher_tmp);
        // cout<<"i:"<<i<<endl;
        // cout<<"cipher1:"<<context.get_context_data(cipher_tmp.parms_id())->chain_index()<<endl;
        evaluator.relinearize_inplace(cipher_tmp, relin_keys);
        evaluator.rescale_to_next_inplace(cipher_tmp);
        // cout<<"cipher1:"<<context.get_context_data(cipher_tmp.parms_id())->chain_index()<<endl;
        // evaluator.rotate_vector_inplace(cipher_tmp, i * rotate_inside_size, galois_keys);
        cipher_tmp_1= rotate_as_16(cipher_tmp, i * rotate_inside_size, scale, encoder, evaluator, decryptor,galois_keys, relin_keys);
        // decrypt_and_print(cipher_tmp, decryptor, encoder);
        // cout<<"cipher1:"<<context.get_context_data(cipher_tmp_1.parms_id())->chain_index()<<endl;
        // cipher_tmp_1.scale()=pow(2.0,45);
        // cout<<"scale:"<<cipher_tmp_1.scale()<<endl;
        cipher_result.push_back(cipher_tmp_1);
    }
    evaluator.add_many(cipher_result, destination);

    
    return destination;
}

seal::Ciphertext pack_change_cipher_to_matrix_A_new_encoding(
    seal::SEALContext context,
    seal::Ciphertext& cipher_matrix, 
    double scale, 
    seal::CKKSEncoder &encoder, 
    seal::Evaluator &evaluator, 
    seal::Decryptor &decryptor,
    seal::GaloisKeys &galois_keys, 
    seal::RelinKeys &relin_keys,
    size_t d)
{
    Matrix<double> u_sigma;
    u_sigma.generate_u_sigma(d, d);
    size_t step=encoder.slot_count()/d/d;
    // cout<<"stepA"<<step<<endl;
    seal::Ciphertext cipher_tmp;//save intermediate variables.
    vector<seal::Ciphertext> cipher_rotate,cipher_result;//save rotate result vector
    seal::Plaintext plain_tmp;
    vector<double> vec_tmp;
    seal::Ciphertext destination;

    //step1. determining the dimensions
    int rotate_size = (2 * d - 1);
    int rotate_outside_size = ceil(sqrt(rotate_size));
    int rotate_inside_size = ceil(double(rotate_size) / double(rotate_outside_size));
    // cout<<rotate_inside_size<<rotate_outside_size<<endl;

    // step2. pre-rotate ciphertext
    

    for (int i = 0; i < rotate_inside_size; i++) {

        evaluator.rotate_vector(cipher_matrix, (-d+d*d+i+1)*step, galois_keys, cipher_tmp);
        // cout<<-d+d*d+i+1<<endl;
        // decrypt_and_print_new_encoding(cipher_tmp, decryptor, encoder);
        cipher_rotate.push_back(cipher_tmp);
    }

    //step3. starting step function 
    for (int i = 0; i < rotate_outside_size; i++) {
        vector<seal::Ciphertext> cipher_vector;
        for (int j = 0; j < rotate_inside_size; j++) {
            vec_tmp = u_sigma.diag_vector(i * rotate_inside_size + j-d+1,d*d);
            // for(size_t s=0;s<vec_tmp.size();s++)
            // {
            //     cout<<vec_tmp[s]<<" ";
            // }
            // cout<<endl;
            // cout<<"vec_tmp.size:"<<vec_tmp.size()<<endl;
            if (std::all_of(vec_tmp.begin(), vec_tmp.end(), [](double num) { return num == 0; }) == 1) {
                continue;
            }
            else {
                std::rotate(vec_tmp.rbegin(), vec_tmp.rbegin() + i * rotate_inside_size, vec_tmp.rend());
                // for(size_t s=0;s<vec_tmp.size();s++)
                // {
                //     cout<<vec_tmp[s]<<" ";
                // }
                // cout<<endl;
                vector<double> result;
                for(size_t s=0;s<vec_tmp.size();s++)
                {
                    for (size_t t=0;t<step;t++){
                        result.push_back(vec_tmp[s]);
                    }
                }
                // 将结果保存回 vec_tmp
                vec_tmp = result;
                encoder.encode(vec_tmp, cipher_rotate[0].scale(), plain_tmp);
                evaluator.mod_switch_to_inplace(plain_tmp, cipher_rotate[0].parms_id());
                evaluator.multiply_plain(cipher_rotate[j], plain_tmp, cipher_tmp);
                
                cipher_vector.push_back(cipher_tmp);
                
            }
        }
        evaluator.add_many(cipher_vector, cipher_tmp);
        // cout<<"i:"<<i<<endl;
        // cout<<"cipher1:"<<context.get_context_data(cipher_tmp.parms_id())->chain_index()<<endl;
        evaluator.relinearize_inplace(cipher_tmp, relin_keys);
        evaluator.rescale_to_next_inplace(cipher_tmp);
        // cout<<"cipher1:"<<context.get_context_data(cipher_tmp.parms_id())->chain_index()<<endl;
        // evaluator.rotate_vector_inplace(cipher_tmp, i * rotate_inside_size, galois_keys);
        evaluator.rotate_vector_inplace(cipher_tmp, (i * rotate_inside_size)*step, galois_keys);
        // decrypt_and_print(cipher_tmp, decryptor, encoder);
        // cout<<"cipher1:"<<context.get_context_data(cipher_tmp_1.parms_id())->chain_index()<<endl;
        // cipher_tmp_1.scale()=pow(2.0,45);
        // cout<<"scale:"<<cipher_tmp_1.scale()<<endl;
        cipher_result.push_back(cipher_tmp);
    }
    evaluator.add_many(cipher_result, destination);

    
    return destination;
}
seal::Ciphertext pack_change_cipher_to_matrix_B(
    seal::SEALContext context,
    seal::Ciphertext& cipher_matrix, 
    double scale, 
    seal::CKKSEncoder &encoder, 
    seal::Evaluator &evaluator, 
    seal::Decryptor &decryptor,
    seal::GaloisKeys &galois_keys, 
    seal::RelinKeys &relin_keys)
{
    // cout<<"start change cipher to matrix B"<<endl;
    // generate u_tau matrix
    size_t d=4;
    Matrix<double> u_tau;
    u_tau.generate_u_tau(d, d);
    seal::Ciphertext cipher_tmp,cipher_tmp_1;//save intermediate variables
    vector<seal::Ciphertext> cipher_rotate, cipher_result;//save rotate result vector
    seal::Plaintext plain_tmp;
    vector<double> vec_tmp;
    seal::Ciphertext destination;

    //start change matrixB
    //step1. determining the dimensions
    int rotate_size = d;
    int rotate_outside_size = ceil(sqrt(rotate_size));
    int rotate_inside_size = ceil(double(rotate_size) / double(rotate_outside_size));
    // step2. pre-rotate ciphertext 
    for (int i = 0; i < rotate_inside_size; i++) {
        cipher_tmp= rotate_as_16(cipher_matrix,i*d, scale, encoder, evaluator, decryptor,galois_keys, relin_keys);
        // decrypt_and_print(cipher_tmp, decryptor, encoder);
        // cipher_tmp.scale()=pow(2.0,45);
        cipher_rotate.push_back(cipher_tmp);
    }
    
    //step3. starting step function 
    for (int i = 0; i < rotate_outside_size; i++) {
        vector<seal::Ciphertext> cipher_vector;
        for (int j = 0; j < rotate_inside_size; j++) {
            if(size_t(i*rotate_inside_size+j)>=d){
                continue;
            }
            vec_tmp = u_tau.diag_vector((i*rotate_inside_size+j)*d,d*d);

            // for(size_t k=0;k<vec_tmp.size();k++){
            //     cout<<vec_tmp[k]<<" ";
            // }
            // cout<<endl;
            if (std::all_of(vec_tmp.begin(), vec_tmp.end(), [](double num) { return num == 0; }) == 1) {
                continue;
            }
            else {
                std::rotate(vec_tmp.rbegin(), vec_tmp.rbegin() + i * rotate_inside_size*d, vec_tmp.rend());
                vector<double> result;
                // 重复 1024 次
                for (int s = 0; s < 1024; s++) {
                    result.insert(result.end(), vec_tmp.begin(), vec_tmp.end());
                }
                // 将结果保存回 vec_tmp
                vec_tmp = result;
                encoder.encode(vec_tmp, cipher_rotate[0].scale(), plain_tmp);
                evaluator.mod_switch_to_inplace(plain_tmp, cipher_rotate[0].parms_id());
                evaluator.multiply_plain(cipher_rotate[j], plain_tmp, cipher_tmp);
                
                cipher_vector.push_back(cipher_tmp);
            }
        }
        evaluator.add_many(cipher_vector, cipher_tmp);
        // cout<<"cipher1:"<<context.get_context_data(cipher_tmp.parms_id())->chain_index()<<endl;
        evaluator.relinearize_inplace(cipher_tmp, relin_keys);
        evaluator.rescale_to_next_inplace(cipher_tmp);
        // cout<<"cipher1:"<<context.get_context_data(cipher_tmp.parms_id())->chain_index()<<endl;
        // evaluator.rotate_vector_inplace(cipher_tmp, i * rotate_inside_size*d, galois_keys);
        
        cipher_tmp_1= rotate_as_16(cipher_tmp, i * rotate_inside_size*d, scale, encoder, evaluator, decryptor,galois_keys, relin_keys);
        // decrypt_and_print(cipher_tmp_1, decryptor, encoder);
        // cout<<"cipher1:"<<context.get_context_data(cipher_tmp_1.parms_id())->chain_index()<<endl;
        // cipher_tmp_1.scale()=pow(2.0,45);
        // cout<<"scale:"<<cipher_tmp_1.scale()<<endl;
        cipher_result.push_back(cipher_tmp_1);
    }
    evaluator.add_many(cipher_result, destination);
    // cipher_matrix=destination;

    return destination;
}
seal::Ciphertext pack_change_cipher_to_matrix_B_new_encoding(
    seal::SEALContext context,
    seal::Ciphertext& cipher_matrix, 
    double scale, 
    seal::CKKSEncoder &encoder, 
    seal::Evaluator &evaluator, 
    seal::Decryptor &decryptor,
    seal::GaloisKeys &galois_keys, 
    seal::RelinKeys &relin_keys,
    size_t d)
{
    Matrix<double> u_tau;
    size_t step=encoder.slot_count()/d/d;
    u_tau.generate_u_tau(d, d);
    seal::Ciphertext cipher_tmp;//save intermediate variables
    vector<seal::Ciphertext> cipher_rotate, cipher_result;//save rotate result vector
    seal::Plaintext plain_tmp;
    vector<double> vec_tmp;
    seal::Ciphertext destination;

    //start change matrixB
    //step1. determining the dimensions
    int rotate_size = d;
    int rotate_outside_size = ceil(sqrt(rotate_size));
    int rotate_inside_size = ceil(double(rotate_size) / double(rotate_outside_size));
    // step2. pre-rotate ciphertext 
    for (int i = 0; i < rotate_inside_size; i++) {
        evaluator.rotate_vector(cipher_matrix, i*d*step, galois_keys, cipher_tmp);
        // decrypt_and_print(cipher_tmp, decryptor, encoder);
        // cipher_tmp.scale()=pow(2.0,45);
        cipher_rotate.push_back(cipher_tmp);
    }
    
    //step3. starting step function 
    for (int i = 0; i < rotate_outside_size; i++) {
        vector<seal::Ciphertext> cipher_vector;
        for (int j = 0; j < rotate_inside_size; j++) {
            if(size_t(i*rotate_inside_size+j)>=d){
                continue;
            }
            vec_tmp = u_tau.diag_vector((i*rotate_inside_size+j)*d,d*d);

            // for(size_t k=0;k<vec_tmp.size();k++){
            //     cout<<vec_tmp[k]<<" ";
            // }
            // cout<<endl;
            if (std::all_of(vec_tmp.begin(), vec_tmp.end(), [](double num) { return num == 0; }) == 1) {
                continue;
            }
            else {
                std::rotate(vec_tmp.rbegin(), vec_tmp.rbegin() + i * rotate_inside_size*d, vec_tmp.rend());
                vector<double> result;
                for(size_t s=0;s<vec_tmp.size();s++)
                {
                    for (size_t t=0;t<step;t++){
                        result.push_back(vec_tmp[s]);
                    }
                }
                // 将结果保存回 vec_tmp
                vec_tmp = result;
                encoder.encode(vec_tmp, cipher_rotate[0].scale(), plain_tmp);
                evaluator.mod_switch_to_inplace(plain_tmp, cipher_rotate[0].parms_id());
                evaluator.multiply_plain(cipher_rotate[j], plain_tmp, cipher_tmp);
                
                cipher_vector.push_back(cipher_tmp);
            }
        }
        evaluator.add_many(cipher_vector, cipher_tmp);
        // cout<<"cipher1:"<<context.get_context_data(cipher_tmp.parms_id())->chain_index()<<endl;
        evaluator.relinearize_inplace(cipher_tmp, relin_keys);
        evaluator.rescale_to_next_inplace(cipher_tmp);
        // cout<<"cipher1:"<<context.get_context_data(cipher_tmp.parms_id())->chain_index()<<endl;
        // evaluator.rotate_vector_inplace(cipher_tmp, i * rotate_inside_size*d, galois_keys);
        
        evaluator.rotate_vector_inplace(cipher_tmp, i * rotate_inside_size*d*step, galois_keys);
        // decrypt_and_print(cipher_tmp_1, decryptor, encoder);
        // cout<<"cipher1:"<<context.get_context_data(cipher_tmp_1.parms_id())->chain_index()<<endl;
        // cipher_tmp_1.scale()=pow(2.0,45);
        // cout<<"scale:"<<cipher_tmp_1.scale()<<endl;
        cipher_result.push_back(cipher_tmp);
    }
    evaluator.add_many(cipher_result, destination);
    // cipher_matrix=destination;

    return destination;
}
seal::Ciphertext rotate_ctA(seal::Ciphertext& cipher_matrix,int i, double scale, 
    seal::CKKSEncoder &encoder, 
    seal::Evaluator &evaluator, 
    seal::Decryptor &decryptor,
    seal::GaloisKeys &galois_keys, 
    seal::RelinKeys &relin_keys)
{
    seal::Ciphertext rotate_result;
    vector<double> vec_tmp(16,0);
    if(i==0){
        vector<double> rotate_all(16384, 1);
        seal::Plaintext plain_tmp;
        encoder.encode(rotate_all, cipher_matrix.scale(), plain_tmp);
        evaluator.mod_switch_to_inplace(plain_tmp, cipher_matrix.parms_id());
        evaluator.multiply_plain(cipher_matrix, plain_tmp, rotate_result);
        evaluator.relinearize_inplace(rotate_result, relin_keys);
        evaluator.rescale_to_next_inplace(rotate_result);



        // rotate_result=cipher_matrix;
        // evaluator.mod_switch_to_next_inplace(rotate_result);

        return rotate_result;
    }
    else{
        size_t d=4;
        seal::Plaintext right_plain, left_plain;
        seal::Ciphertext right_cipher,right_cipher_tmp, left_cipher,left_cipher_tmp;
        vector<double> rotate_right(d * d), rotate_left(d * d);
        int w=d-i;
        for (size_t j = 0; j < d; j++) {
            size_t k = 0;
            while (k < d) {
                if (k >= d-w) {
                    rotate_right[j * d + k] = 1;
                }
                else {
                    rotate_left[j * d + k] = 1;
                }
                k++;
            }
        }
        // for(int k=0;k<16;k++)
        // {
        //     cout<<rotate_right[k]<<" ";
        // }
        // cout<<endl;
        vector<double> result;
        // 重复 1024 次
        for (int s = 0; s < 1024; s++) {
            result.insert(result.end(), rotate_left.begin(), rotate_left.end());
        }
        rotate_left = result;
        encoder.encode(rotate_left,cipher_matrix.scale(),left_plain);
        evaluator.mod_switch_to_inplace(left_plain, cipher_matrix.parms_id());
        evaluator.multiply_plain(cipher_matrix, left_plain, left_cipher_tmp);
        evaluator.rotate_vector(left_cipher_tmp, i-4, galois_keys, left_cipher);
        // decrypt_and_print(left_cipher, decryptor, encoder);

        vector<double> result1;
        // 重复 1024 次
        for (int s = 0; s < 1024; s++) {
            result1.insert(result1.end(), rotate_right.begin(),rotate_right.end());
        }
        rotate_right = result1;
        encoder.encode(rotate_right, cipher_matrix.scale(), right_plain);
        evaluator.mod_switch_to_inplace(right_plain, cipher_matrix.parms_id());
        evaluator.multiply_plain(cipher_matrix, right_plain, right_cipher_tmp);
        evaluator.rotate_vector(right_cipher_tmp, i, galois_keys, right_cipher);
        // decrypt_and_print(right_cipher, decryptor, encoder);
        evaluator.add(left_cipher,right_cipher,rotate_result);
        evaluator.relinearize_inplace(rotate_result, relin_keys);
        evaluator.rescale_to_next_inplace(rotate_result);
        return rotate_result;
    }
}
seal::Ciphertext pack_encrypted_matrix_mul(
    seal::SEALContext context,
    seal::Ciphertext& cipher_A, 
    seal::Ciphertext& cipher_B, 
    double scale, 
    seal::CKKSEncoder &encoder, 
    seal::Evaluator &evaluator, 
    seal::Decryptor &decryptor,
    seal::GaloisKeys &galois_keys, 
    seal::RelinKeys &relin_keys)
{
    size_t d=4;
    vector<seal::Ciphertext>  cipher_vector;
    seal::Ciphertext cipher_tmpA,cipher_tmpB,cipher_tmp,encrypted_AB;//save intermediate variables
    seal::Plaintext plain_tmp;
    vector<double> vec_tmp;
    for (size_t i = 0; i < d; i++) {
        cipher_tmpA=rotate_ctA(cipher_A, i, scale, encoder, evaluator, decryptor,galois_keys, relin_keys);
        cipher_tmpB=rotate_as_16(cipher_B,d*i, scale, encoder, evaluator, decryptor,galois_keys, relin_keys);
        evaluator.mod_switch_to_inplace(cipher_tmpB, cipher_tmpA.parms_id());
        evaluator.multiply(cipher_tmpA, cipher_tmpB, cipher_tmp);
        cipher_vector.push_back(cipher_tmp);
    }
    evaluator.add_many(cipher_vector, encrypted_AB);
    evaluator.relinearize_inplace(encrypted_AB, relin_keys);
    evaluator.rescale_to_next_inplace(encrypted_AB);


    return encrypted_AB;
}

seal::Ciphertext rotate_ctA_new_encoding(seal::Ciphertext& cipher_matrix,int i, double scale, 
    seal::CKKSEncoder &encoder, 
    seal::Evaluator &evaluator, 
    seal::Decryptor &decryptor,
    seal::GaloisKeys &galois_keys, 
    seal::RelinKeys &relin_keys,
    size_t d)
{
    seal::Ciphertext rotate_result;
    size_t step=encoder.slot_count()/d/d;
    // vector<double> vec_tmp(16,0);
    if(i==0){
        vector<double> rotate_all(16384, 1);
        seal::Plaintext plain_tmp;
        encoder.encode(rotate_all, cipher_matrix.scale(), plain_tmp);
        evaluator.mod_switch_to_inplace(plain_tmp, cipher_matrix.parms_id());
        evaluator.multiply_plain(cipher_matrix, plain_tmp, rotate_result);
        evaluator.relinearize_inplace(rotate_result, relin_keys);
        evaluator.rescale_to_next_inplace(rotate_result);



        // rotate_result=cipher_matrix;
        // evaluator.mod_switch_to_next_inplace(rotate_result);

        return rotate_result;
    }
    else{
        seal::Plaintext right_plain, left_plain;
        seal::Ciphertext right_cipher,right_cipher_tmp, left_cipher,left_cipher_tmp;
        vector<double> rotate_right(d * d), rotate_left(d * d);
        int w=d-i;
        for (size_t j = 0; j < d; j++) {
            size_t k = 0;
            while (k < d) {
                if (k >= d-w) {
                    rotate_right[j * d + k] = 1;
                }
                else {
                    rotate_left[j * d + k] = 1;
                }
                k++;
            }
        }
        // for(size_t k=0;k<rotate_right.size();k++)
        // {
        //     cout<<rotate_right[k]<<" ";
        // }
        // cout<<endl;
        // cout<<"left"<<endl;
        // for(size_t k=0;k<rotate_left.size();k++)
        // {
        //     cout<<rotate_left[k]<<" ";
        // }
        // cout<<endl;

        vector<double> result;
        for(size_t s=0;s<rotate_left.size();s++)
        {
            for (size_t t=0;t<step;t++){
                result.push_back(rotate_left[s]);
            }
        }

        // 将结果保存回 
        rotate_left = result;

        encoder.encode(rotate_left,cipher_matrix.scale(),left_plain);
        evaluator.mod_switch_to_inplace(left_plain, cipher_matrix.parms_id());
        evaluator.multiply_plain(cipher_matrix, left_plain, left_cipher_tmp);
        evaluator.rotate_vector(left_cipher_tmp, (i-d)*step, galois_keys, left_cipher);
        // decrypt_and_print(left_cipher, decryptor, encoder);

        vector<double> result1;
        for(size_t s=0;s<rotate_right.size();s++)
        {
            for (size_t t=0;t<step;t++){
                result1.push_back(rotate_right[s]);
            }
        }
        // 将结果保存回 
        rotate_right = result1;
        encoder.encode(rotate_right, cipher_matrix.scale(), right_plain);
        evaluator.mod_switch_to_inplace(right_plain, cipher_matrix.parms_id());
        evaluator.multiply_plain(cipher_matrix, right_plain, right_cipher_tmp);
        evaluator.rotate_vector(right_cipher_tmp, step*i, galois_keys, right_cipher);
        // decrypt_and_print(right_cipher, decryptor, encoder);
        evaluator.add(left_cipher,right_cipher,rotate_result);
        evaluator.relinearize_inplace(rotate_result, relin_keys);
        evaluator.rescale_to_next_inplace(rotate_result);
        return rotate_result;
    }
}
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
    size_t d)
{
    size_t step=encoder.slot_count()/d/d;
    vector<seal::Ciphertext>  cipher_vector;
    seal::Ciphertext cipher_tmpA,cipher_tmpB,cipher_tmp,encrypted_AB;//save intermediate variables
    seal::Plaintext plain_tmp;
    vector<double> vec_tmp;
    for (size_t i = 0; i < d; i++) {
        cipher_tmpA=rotate_ctA_new_encoding(cipher_A, i, scale, encoder, evaluator, decryptor,galois_keys, relin_keys,d);
        // decrypt_and_print_new_encoding(cipher_tmpA, decryptor, encoder,d);
        evaluator.rotate_vector(cipher_B,i*d*step,galois_keys,cipher_tmpB);
        evaluator.mod_switch_to_inplace(cipher_tmpB, cipher_tmpA.parms_id());
        evaluator.multiply(cipher_tmpA, cipher_tmpB, cipher_tmp);
        cipher_vector.push_back(cipher_tmp);
    }
    evaluator.add_many(cipher_vector, encrypted_AB);
    evaluator.relinearize_inplace(encrypted_AB, relin_keys);
    evaluator.rescale_to_next_inplace(encrypted_AB);


    return encrypted_AB;
}

