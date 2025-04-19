#include<iostream>
#include <seal/seal.h>
#include "../src/jiang_enc_cipher.h"
#include "../src/block.h"
#include "../src/matrix.h"
#include "../src/Bootstrapper.h"
#include <cstdlib>
#include <chrono>
#include <fstream>

using namespace std;
using namespace seal;
void run_jiang_mat_mul(size_t n, size_t m, size_t p, size_t poly_modulus_degree);
void test_jiang_mat_mul(SEALContext context,size_t n, size_t m, size_t p, 
                        CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                        RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);
void run_dia_block_inversion(size_t n, size_t m, size_t p, size_t poly_modulus_degree);
void test_dia_block_inversion(SEALContext context,size_t n, size_t m, size_t p, 
                        CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                        RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);
void run_real_block_inversion(size_t n, size_t m, size_t p, size_t poly_modulus_degree);
void test_real_block_inversion(SEALContext context,size_t n, size_t m, size_t p, 
                        CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                        RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);
void run_pack_inversion(size_t n, size_t m, size_t p, size_t poly_modulus_degree);
void test_pack_inversion(SEALContext context,size_t n, size_t m, size_t p, 
                        CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                        RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);
void run_dia_pack_new_encoding_inversion(size_t n, size_t m, size_t p, size_t poly_modulus_degree);
void test_dia_pack_new_encoding_inversion(SEALContext context,size_t n, size_t m, size_t p, 
                        CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                        RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);
void run_real_pack_new_encoding_inversion(size_t n, size_t m, size_t p, size_t poly_modulus_degree);
void test_real_pack_new_encoding_inversion(SEALContext context,size_t n, size_t m, size_t p, 
                        CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                        RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);
void run_triangle_inversion(size_t n, size_t m, size_t p,size_t poly_modulus_degree);
void test_triangle_inversion(SEALContext context, size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);
void run_tae_inversion(size_t n, size_t m, size_t p,size_t poly_modulus_degree);
void test_tae_inversion(SEALContext context, size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);
void run_dia_block_inversion_512(size_t n, size_t m, size_t p, size_t poly_modulus_degree);
void test_dia_block_inversion_512(SEALContext context,size_t n, size_t m, size_t p, 
                        CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                        RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);
void run_real_block_inversion_512(size_t n, size_t m, size_t p, size_t poly_modulus_degree);
void test_real_block_inversion_512(SEALContext context,size_t n, size_t m, size_t p, 
                        CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, 
                        RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor);
int main()
{
    srand(123);
    while(true)
    {
        cout << "\nExamples:" << endl
             << endl;

        cout << "  1.  triangle inv(4,4)" << endl;
        cout << "  2.  dia: pack new encoding(d,d)" << endl;
        cout << "  3.  real: pack new encoding(d,d)" << endl;
        cout << "  4.  tae inv(d,d)" << endl;
        cout << "  5.  dia:block inversion(256,256)" << endl;
        cout << "  6.  dia:block inversion(512,512)" << endl;
        cout << "  7-10. test"<<endl;
        
        // cout << "  7.  native inversion(128,128)" << endl;                
        // cout << "  8.  pack inversion with row encoding(4,4)" << endl;
        // cout << "  9.  real:block inversion(256,256)" << endl;
        // cout << "  10. real:block inversion(512,512)" << endl;
        cout << "  0. Exit" << endl;
        
        cout << "\nTotal memory allocated from the current memory pool: "
             << (MemoryManager::GetPool().alloc_byte_count() >> 20) << " MB" << endl;

        int selection = 0;
        cout << endl
             << "Run example: ";
        if (!(cin >> selection))
        {
            cout << "Invalid option." << endl;
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            continue;
        }
        switch (selection)
        {
            case 7:
            //     run_jiang_mat_mul(64,64,64,8192);
                run_jiang_mat_mul(128,128,128,32768);
                break;
            case 5:
                run_dia_block_inversion(256,256,256,32768);
                break;
            case 8:
                run_pack_inversion(4,4,4,32768);
                break;
            case 2:
                int d;
                cout<<"请输入维数：";
                cin>>d;
                run_dia_pack_new_encoding_inversion(d,d,d,32768);
                break;
            case 3:
                int d_real;
                cout<<"请输入维数：";
                cin>>d_real;
                run_real_pack_new_encoding_inversion(d_real,d_real,d_real,32768);
                break;
            case 1:
                size_t N;
                cout<<"请输入参数N:";
                cin>>N;
                run_triangle_inversion(4,4,4,N);
                break;
            case 4:
                int d1;
                cout<<"请输入维数：";
                cin>>d1;
                run_tae_inversion(d1,d1,d1,32768);
                break;
            case 6:
                run_dia_block_inversion_512(512,512,512,32768);
                break;
            case 9:
                run_real_block_inversion(256,256,256,32768);
                break;
            case 10:
                run_real_block_inversion_512(512,512,512,32768);
                break;
            case 0:
                return 0;
            default:
                cout<<"Invalid option." << endl;
        }
    }


    return 0;
}
void run_jiang_mat_mul(size_t n, size_t m, size_t p,size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
//     parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30,30,30, 60}));
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 40, 40, 40, 40, 40, 40, 40, 40, 60}));
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 45, 45, 45, 45, 45, 45, 45, 45, 60}));
    double scale = pow(2.0, 45);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    // print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);
    cout<<"---------"<<endl;
    test_jiang_mat_mul(context,n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}
void test_jiang_mat_mul(SEALContext context, size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
     /*初始化矩阵*/
     Matrix<double> A(n, m);
     // Matrix<double> B(m, p);
     size_t d=sqrt(encoder.slot_count());
     random_matrix_generator(n,m,A);

     Matrix<double> I;
     eye_matrix_generator(n,m,I);

     
     Matrix<double> C;

     C = A.multiply(A);

     Matrix<double> D;

     D = C.add(A);

     Matrix<double> E;
     E = C.multiply(D);
     
     Matrix<double> F;
     F = D.add(E);

     auto start_time = chrono::high_resolution_clock::now();

     cipher_matrix_jiang cipher_A0(A,scale,encoder,encryptor);//A
     seal::Ciphertext cipher_tmps0;
     cipher_tmps0=cipher_A0.get_cipher_matrix();
    
     seal::Plaintext plain0,plain_I;
     encoder.encode(-0.1,cipher_tmps0.scale(),plain0);
     evaluator.mod_switch_to_inplace(plain0,cipher_tmps0.parms_id());
     evaluator.multiply_plain_inplace(cipher_tmps0,plain0);


     vector<double> flatten_data=I.flatten_matrix_to_rows_vector();
     encoder.encode(flatten_data,cipher_tmps0.scale(),plain_I);
     evaluator.add_plain_inplace(cipher_tmps0,plain_I);

     evaluator.rescale_to_next_inplace(cipher_tmps0);
     cipher_matrix_jiang cipher_A(cipher_tmps0,n,m,d);




     cipher_matrix_jiang cipher_B;//A
     cipher_matrix_jiang cipher_A_square;//A^2

     cipher_matrix_jiang cipher_A_fourth;//A^2(A+A^2)

     seal::Ciphertext cipher_A_ciphertext;
     cipher_A_ciphertext=cipher_A.get_cipher_matrix();

     cipher_B=cipher_A;
     cipher_A.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
     cipher_B.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
     encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A,cipher_B,cipher_A_square);
     
     seal::Ciphertext cipher_A_ciphertext_2,cipher_tmps;
     cipher_A_ciphertext_2=cipher_A_square.get_cipher_matrix();
     match_levels(context,cipher_A_ciphertext, cipher_A_ciphertext_2, evaluator);
     evaluator.add(cipher_A_ciphertext,cipher_A_ciphertext_2,cipher_tmps);
     cipher_matrix_jiang cipher_C(cipher_tmps,n,m,d);//A^2+A

     seal::Ciphertext cipher_C_ciphertext;
     cipher_C_ciphertext=cipher_C.get_cipher_matrix();

     cipher_A_square.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
     cipher_C.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
     encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A_square,cipher_C,cipher_A_fourth);
     
     seal::Ciphertext cipher_A_ciphertext_4,cipher_tmps1;
     cipher_A_ciphertext_4=cipher_A_fourth.get_cipher_matrix();
     match_levels(context,cipher_C_ciphertext, cipher_A_ciphertext_4, evaluator);
     evaluator.add(cipher_C_ciphertext,cipher_A_ciphertext_4,cipher_tmps1);

     seal::Plaintext plain1;
     encoder.encode(0.1,cipher_tmps1.scale(),plain1);
     evaluator.mod_switch_to_inplace(plain1,cipher_tmps1.parms_id());
     evaluator.multiply_plain_inplace(cipher_tmps1,plain1);

    //  seal::Plaintext plain_I_2;
    //  vector<double> flatten_data2=I.flatten_matrix_to_rows_vector();
    //  encoder.encode(flatten_data2,cipher_tmps1.scale(),plain_I_2);
    //  evaluator.add_plain_inplace(cipher_tmps1,plain_I_2);
     cipher_matrix_jiang result(cipher_tmps1,n,m,d);//A+A^2+A^3+A^4


     // std::cout << "Noise budget in encrypted: " << decryptor.invariant_noise_budget(cipher_tmps) << " bits" << std::endl;
     // // cout<<"cipher1:"<<context.get_context_data(cipher_A.get_cipher_matrix().parms_id())->chain_index()<<endl;
     // // cout<<log2(cipher_A.get_cipher_matrix().scale())<<"bits"<<endl; 




     // decrypte result

     Matrix<double> computed_C;
     result.dec_matrix_cipher(computed_C,encoder,decryptor);
     auto end_time = chrono::high_resolution_clock::now();
     auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
     cout << endl
          << "----------------- Jiang matrix multiplication --------------------" << endl;
     cout << endl
          << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
     cout << "Actual Multiply:"
          << " (" << d << "x" << d << ")X(" << d << "x" << d << ")" << endl;
     cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

     cout << "the computed C = A^-1: " << endl;
    //  store_matrix(n,m,computed_C);
     computed_C.print();
         // // 将矩阵存储到文件
    //  std::ofstream outFile("matrix_result_128.txt");
    //  for (size_t i = 0; i < m; i++) {
    //      for (size_t j = 0; j < n; j++) {
    //          outFile << computed_C.get(j, i) << (j < n - 1 ? " " : "");
    //      }
    //      outFile << std::endl;
    //  }

    //  outFile.close();

     cout<<"endl"<<"---------------------------------------------------------"<<endl;

     cout << "the true C =  A^2+A: " << endl;
     F.print();
     computed_C.scaling_inplace(-1.);
     computed_C = computed_C.add(C);
     computed_C.is_zero(-2);
     cout << "------------------------------------------------------------------" << endl;
     cout << endl <<endl;
}

void run_dia_block_inversion(size_t n, size_t m, size_t p,size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
//     parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30,30,30, 60}));
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 40, 40, 40, 40, 40, 40, 40, 40, 60}));
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 45, 45, 45, 45, 45, 45, 45, 45, 60}));
    double scale = pow(2.0, 45);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    // print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);
    cout<<"---------"<<endl;
    test_dia_block_inversion(context,n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}
void test_dia_block_inversion(SEALContext context,size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    // Matrix<double>A(n,m);
    // random_block_matrix_generator(n, m, A);
    // Matrix<double> C;//A^2
    // C = A.multiply(A);

    // Matrix<double>D;//A^2+A
    // D = C.add(A);
    size_t d=sqrt(encoder.slot_count());

    // Matrix<double>E;//A^2+A
    // E = C.multiply(D);

    // Matrix<double>F;//A^2+A+A^2+A
    // F = E.add(D);



    // Matrix<double>A_true_0(d,d);
    // Matrix<double>A_true_1(d,d);
    // Matrix<double>A_true_2(d,d);
    // Matrix<double>A_true_3(d,d);
    // F.split(d,A_true_0,A_true_1,A_true_2,A_true_3);
    Matrix<double>A(n,m);
    std::string filename = "../build/dia_input/dia_input_c_(256,256).txt";
    A.readFromFile(filename);
    cout<<"   wssssssssss"<<endl;
    // cout<<A.get(0,1);

    A.print(); // 打印矩阵以验证




    Matrix<double>A_0(d,d);
    Matrix<double>A_1(d,d);
    Matrix<double>A_2(d,d);
    Matrix<double>A_3(d,d);

    A.split(d,A_0, A_1,A_2,A_3);


    auto start_time = chrono::high_resolution_clock::now();



    cipher_matrix_jiang cipher_A0(A_0,scale,encoder,encryptor);//A_0
    cipher_matrix_jiang cipher_A1(A_1,scale,encoder,encryptor);//A_1
    cipher_matrix_jiang cipher_A2(A_2,scale,encoder,encryptor);//A_2
    cipher_matrix_jiang cipher_A3(A_3,scale,encoder,encryptor);//A_3

    cipher_matrix_jiang cipher_B0;
    cipher_matrix_jiang cipher_B1;
    cipher_matrix_jiang cipher_B2;
    cipher_matrix_jiang cipher_B3;
    cipher_B0=cipher_A0;
    cipher_B1=cipher_A1;
    cipher_B2=cipher_A2;
    cipher_B3=cipher_A3;


    seal::Ciphertext cipher_A0_ciphertext,cipher_A1_ciphertext,cipher_A2_ciphertext,cipher_A3_ciphertext;
    cipher_A0_ciphertext=cipher_A0.get_cipher_matrix();
    cipher_A1_ciphertext=cipher_A1.get_cipher_matrix();
    cipher_A2_ciphertext=cipher_A2.get_cipher_matrix();
    cipher_A3_ciphertext=cipher_A3.get_cipher_matrix();

    cipher_A0.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B0.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A1.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B1.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A2.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B2.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A3.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B3.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    

    //compute A^2
    cipher_matrix_jiang cipher_A0B0;
    cipher_matrix_jiang cipher_A1B2;
    cipher_matrix_jiang cipher_A0B1;
    cipher_matrix_jiang cipher_A1B3;
    cipher_matrix_jiang cipher_A2B0;
    cipher_matrix_jiang cipher_A3B2;
    cipher_matrix_jiang cipher_A2B1;
    cipher_matrix_jiang cipher_A3B3;
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A0,cipher_B0,cipher_A0B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1,cipher_B2,cipher_A1B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A0,cipher_B1,cipher_A0B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1,cipher_B3,cipher_A1B3);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A2,cipher_B0,cipher_A2B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A3,cipher_B2,cipher_A3B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A2,cipher_B1,cipher_A2B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A3,cipher_B3,cipher_A3B3);
    
    seal::Ciphertext cipher_A0B0_ciphertext,cipher_A1B2_ciphertext,cipher_A0B1_ciphertext,cipher_A1B3_ciphertext,cipher_A2B0_ciphertext,cipher_A3B2_ciphertext,cipher_A2B1_ciphertext,cipher_A3B3_ciphertext;
    cipher_A0B0_ciphertext=cipher_A0B0.get_cipher_matrix();
    cipher_A1B2_ciphertext=cipher_A1B2.get_cipher_matrix();
    cipher_A0B1_ciphertext=cipher_A0B1.get_cipher_matrix();
    cipher_A1B3_ciphertext=cipher_A1B3.get_cipher_matrix();
    cipher_A2B0_ciphertext=cipher_A2B0.get_cipher_matrix();
    cipher_A3B2_ciphertext=cipher_A3B2.get_cipher_matrix();
    cipher_A2B1_ciphertext=cipher_A2B1.get_cipher_matrix();
    cipher_A3B3_ciphertext=cipher_A3B3.get_cipher_matrix();

    //A^2 
    seal::Ciphertext cipher_temp_0,cipher_temp_1,cipher_temp_2,cipher_temp_3;
    evaluator.add(cipher_A0B0_ciphertext,cipher_A1B2_ciphertext,cipher_temp_0);
    evaluator.add(cipher_A0B1_ciphertext,cipher_A1B3_ciphertext,cipher_temp_1);
    evaluator.add(cipher_A2B0_ciphertext,cipher_A3B2_ciphertext,cipher_temp_2);
    evaluator.add(cipher_A2B1_ciphertext,cipher_A3B3_ciphertext,cipher_temp_3);
    cipher_matrix_jiang A_square_0(cipher_temp_0,d,d,d);
    cipher_matrix_jiang A_square_1(cipher_temp_1,d,d,d);
    cipher_matrix_jiang A_square_2(cipher_temp_2,d,d,d);
    cipher_matrix_jiang A_square_3(cipher_temp_3,d,d,d);
    
    //compute A^2+A
    seal::Ciphertext cipher_temp_0_A0,cipher_temp_1_A1,cipher_temp_2_A2,cipher_temp_3_A3;
    match_levels(context,cipher_A0_ciphertext,cipher_temp_0,evaluator);
    evaluator.add(cipher_temp_0,cipher_A0_ciphertext,cipher_temp_0_A0);

    match_levels(context,cipher_A1_ciphertext,cipher_temp_1,evaluator);
    evaluator.add(cipher_temp_1,cipher_A1_ciphertext,cipher_temp_1_A1);

    match_levels(context,cipher_A2_ciphertext,cipher_temp_2,evaluator);
    evaluator.add(cipher_temp_2,cipher_A2_ciphertext,cipher_temp_2_A2);

    match_levels(context,cipher_A3_ciphertext,cipher_temp_3,evaluator);
    evaluator.add(cipher_temp_3,cipher_A3_ciphertext,cipher_temp_3_A3);

    cipher_matrix_jiang A_temp_0(cipher_temp_0_A0,d,d,d);
    cipher_matrix_jiang A_temp_1(cipher_temp_1_A1,d,d,d);
    cipher_matrix_jiang A_temp_2(cipher_temp_2_A2,d,d,d);
    cipher_matrix_jiang A_temp_3(cipher_temp_3_A3,d,d,d);
    //compute A^2(A^2+A)

    A_square_0.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    A_temp_0.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    A_square_1.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    A_temp_1.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    A_square_2.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    A_temp_2.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    A_square_3.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    A_temp_3.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);

    cipher_matrix_jiang cipher_A0B0_square;
    cipher_matrix_jiang cipher_A1B2_square;
    cipher_matrix_jiang cipher_A0B1_square;
    cipher_matrix_jiang cipher_A1B3_square;
    cipher_matrix_jiang cipher_A2B0_square;
    cipher_matrix_jiang cipher_A3B2_square;
    cipher_matrix_jiang cipher_A2B1_square;
    cipher_matrix_jiang cipher_A3B3_square;
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,A_square_0,A_temp_0,cipher_A0B0_square);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,A_square_1,A_temp_2,cipher_A1B2_square);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,A_square_0,A_temp_1,cipher_A0B1_square);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,A_square_1,A_temp_3,cipher_A1B3_square);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,A_square_2,A_temp_0,cipher_A2B0_square);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,A_square_3,A_temp_2,cipher_A3B2_square);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,A_square_2,A_temp_1,cipher_A2B1_square);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,A_square_3,A_temp_3,cipher_A3B3_square);
    
    
    seal::Ciphertext cipher_A0B0_ciphertext_square,cipher_A1B2_ciphertext_square,cipher_A0B1_ciphertext_square,cipher_A1B3_ciphertext_square,cipher_A2B0_ciphertext_square,cipher_A3B2_ciphertext_square,cipher_A2B1_ciphertext_square,cipher_A3B3_ciphertext_square;
    cipher_A0B0_ciphertext_square=cipher_A0B0_square.get_cipher_matrix();
    cipher_A1B2_ciphertext_square=cipher_A1B2_square.get_cipher_matrix();
    cipher_A0B1_ciphertext_square=cipher_A0B1_square.get_cipher_matrix();
    cipher_A1B3_ciphertext_square=cipher_A1B3_square.get_cipher_matrix();
    cipher_A2B0_ciphertext_square=cipher_A2B0_square.get_cipher_matrix();
    cipher_A3B2_ciphertext_square=cipher_A3B2_square.get_cipher_matrix();
    cipher_A2B1_ciphertext_square=cipher_A2B1_square.get_cipher_matrix();
    cipher_A3B3_ciphertext_square=cipher_A3B3_square.get_cipher_matrix();

    
    seal::Ciphertext cipher_temp_0_square,cipher_temp_1_square,cipher_temp_2_square,cipher_temp_3_square;
    evaluator.add(cipher_A0B0_ciphertext_square,cipher_A1B2_ciphertext_square,cipher_temp_0_square);
    evaluator.add(cipher_A0B1_ciphertext_square,cipher_A1B3_ciphertext_square,cipher_temp_1_square);
    evaluator.add(cipher_A2B0_ciphertext_square,cipher_A3B2_ciphertext_square,cipher_temp_2_square);
    evaluator.add(cipher_A2B1_ciphertext_square,cipher_A3B3_ciphertext_square,cipher_temp_3_square);


    //compute (A^2+A)  with (A^3+A^4)
    seal::Ciphertext result0,result1,result2,result3;
    match_levels(context,cipher_temp_0_A0,cipher_temp_0_square,evaluator);
    evaluator.add(cipher_temp_0_A0,cipher_temp_0_square,result0);
    match_levels(context,cipher_temp_1_A1,cipher_temp_1_square,evaluator);
    evaluator.add(cipher_temp_1_A1,cipher_temp_1_square,result1);
    match_levels(context,cipher_temp_2_A2,cipher_temp_2_square,evaluator);
    evaluator.add(cipher_temp_2_A2,cipher_temp_2_square,result2);
    match_levels(context,cipher_temp_3_A3,cipher_temp_3_square,evaluator);
    evaluator.add(cipher_temp_3_A3,cipher_temp_3_square,result3);

    cipher_matrix_jiang result_0(result0,d,d,d);
    cipher_matrix_jiang result_1(result1,d,d,d);
    cipher_matrix_jiang result_2(result2,d,d,d);
    cipher_matrix_jiang result_3(result3,d,d,d);



    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    Matrix<double> computed0,computed1,computed2,computed3;
    result_0.dec_matrix_cipher(computed0,encoder,decryptor);
    result_1.dec_matrix_cipher(computed1,encoder,decryptor);
    result_2.dec_matrix_cipher(computed2,encoder,decryptor);
    result_3.dec_matrix_cipher(computed3,encoder,decryptor);

    std::ofstream outFile0("matrix_result_256_0.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile0 << computed0.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile0 << std::endl;
    }

    outFile0.close();
    std::ofstream outFile1("matrix_result_256_1.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile1 << computed1.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile1 << std::endl;
    }

    outFile1.close();
    std::ofstream outFile2("matrix_result_256_2.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile2 << computed2.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile2 << std::endl;
    }

    outFile2.close();
    std::ofstream outFile3("matrix_result_256_3.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile3 << computed3.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile3 << std::endl;
    }

    outFile3.close();

    cout<<"endl"<<"---------------------------------------------------------"<<endl;
    // cout << "the true C =  A^2+A: " << endl;
    // A_true_0.print();

}

void run_real_block_inversion(size_t n, size_t m, size_t p,size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 60}));
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 45, 45, 45, 45, 45, 45, 45, 45, 60}));
    double scale = pow(2.0, 45);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    // print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);
    cout<<"---------"<<endl;
    test_real_block_inversion(context,n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}
void test_real_block_inversion(SEALContext context,size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    // Matrix<double>A(n,m);
    // random_block_matrix_generator(n, m, A);
    // Matrix<double> C;//A^2
    // C = A.multiply(A);

    // Matrix<double>D;//A^2+A
    // D = C.add(A);
    size_t d=sqrt(encoder.slot_count());

    // Matrix<double>E;//A^2+A
    // E = C.multiply(D);

    // Matrix<double>F;//A^2+A+A^2+A
    // F = E.add(D);



    // Matrix<double>A_true_0(d,d);
    // Matrix<double>A_true_1(d,d);
    // Matrix<double>A_true_2(d,d);
    // Matrix<double>A_true_3(d,d);
    // F.split(d,A_true_0,A_true_1,A_true_2,A_true_3);
    Matrix<double>A(n,m);
    std::string filename = "../build/real_input/real_input_c_(256,256).txt";
    A.readFromFile(filename);
    cout<<"   wssssssssss"<<endl;
    // cout<<A.get(0,1);

    A.print(); // 打印矩阵以验证




    Matrix<double>A_0(d,d);
    Matrix<double>A_1(d,d);
    Matrix<double>A_2(d,d);
    Matrix<double>A_3(d,d);

    A.split(d,A_0, A_1,A_2,A_3);


    auto start_time = chrono::high_resolution_clock::now();



    cipher_matrix_jiang cipher_A0(A_0,scale,encoder,encryptor);//A_0
    cipher_matrix_jiang cipher_A1(A_1,scale,encoder,encryptor);//A_1
    cipher_matrix_jiang cipher_A2(A_2,scale,encoder,encryptor);//A_2
    cipher_matrix_jiang cipher_A3(A_3,scale,encoder,encryptor);//A_3

    cipher_matrix_jiang cipher_B0;
    cipher_matrix_jiang cipher_B1;
    cipher_matrix_jiang cipher_B2;
    cipher_matrix_jiang cipher_B3;
    cipher_B0=cipher_A0;
    cipher_B1=cipher_A1;
    cipher_B2=cipher_A2;
    cipher_B3=cipher_A3;
    // A^2
    vector<Ciphertext> left_degree_2=block_256_LeftSquare(context,encoder,scale,evaluator,relin_keys,gal_keys,cipher_A0,cipher_A1,cipher_A2,cipher_A3,cipher_B0,cipher_B1,cipher_B2,cipher_B3);
    //A+A*A
    vector<Ciphertext> right_degree_2=block_256_LeftMulRight_AddRight(context,encoder,scale,evaluator,relin_keys,gal_keys,cipher_A0,cipher_A1,cipher_A2,cipher_A3,cipher_B0,cipher_B1,cipher_B2,cipher_B3);
    cipher_matrix_jiang left_degree_2_0(left_degree_2[0],d,d,d);
    cipher_matrix_jiang left_degree_2_1(left_degree_2[1],d,d,d);
    cipher_matrix_jiang left_degree_2_2(left_degree_2[2],d,d,d);
    cipher_matrix_jiang left_degree_2_3(left_degree_2[3],d,d,d);
    cipher_matrix_jiang right_degree_2_0(right_degree_2[0],d,d,d);
    cipher_matrix_jiang right_degree_2_1(right_degree_2[1],d,d,d);
    cipher_matrix_jiang right_degree_2_2(right_degree_2[2],d,d,d);
    cipher_matrix_jiang right_degree_2_3(right_degree_2[3],d,d,d);


    // A^4
    vector<Ciphertext> left_degree_4=block_256_LeftSquare(context,encoder,scale,evaluator,relin_keys,gal_keys,left_degree_2_0,left_degree_2_1,left_degree_2_2,left_degree_2_3,left_degree_2_0,left_degree_2_1,left_degree_2_2,left_degree_2_3);
    //A+...+A^4
    vector<Ciphertext> right_degree_4=block_256_LeftMulRight_AddRight(context,encoder,scale,evaluator,relin_keys,gal_keys,left_degree_2_0,left_degree_2_1,left_degree_2_2,left_degree_2_3,right_degree_2_0,right_degree_2_1,right_degree_2_2,right_degree_2_3);
    cipher_matrix_jiang left_degree_4_0(left_degree_4[0],d,d,d);
    cipher_matrix_jiang left_degree_4_1(left_degree_4[1],d,d,d);
    cipher_matrix_jiang left_degree_4_2(left_degree_4[2],d,d,d);
    cipher_matrix_jiang left_degree_4_3(left_degree_4[3],d,d,d);
    cipher_matrix_jiang right_degree_4_0(right_degree_4[0],d,d,d);
    cipher_matrix_jiang right_degree_4_1(right_degree_4[1],d,d,d);
    cipher_matrix_jiang right_degree_4_2(right_degree_4[2],d,d,d);
    cipher_matrix_jiang right_degree_4_3(right_degree_4[3],d,d,d);
    vector<Ciphertext> left_degree_8=block_256_LeftSquare(context,encoder,scale,evaluator,relin_keys,gal_keys,left_degree_4_0,left_degree_4_1,left_degree_4_2,left_degree_4_3,left_degree_4_0,left_degree_4_1,left_degree_4_2,left_degree_4_3);
    //A+...+A^8
    vector<Ciphertext> right_degree_8=block_256_LeftMulRight_AddRight(context,encoder,scale,evaluator,relin_keys,gal_keys,left_degree_4_0,left_degree_4_1,left_degree_4_2,left_degree_4_3,right_degree_4_0,right_degree_4_1,right_degree_4_2,right_degree_4_3);
    cipher_matrix_jiang left_degree_8_0(left_degree_8[0],d,d,d);
    cipher_matrix_jiang left_degree_8_1(left_degree_8[1],d,d,d);
    cipher_matrix_jiang left_degree_8_2(left_degree_8[2],d,d,d);
    cipher_matrix_jiang left_degree_8_3(left_degree_8[3],d,d,d);
    cipher_matrix_jiang right_degree_8_0(right_degree_8[0],d,d,d);
    cipher_matrix_jiang right_degree_8_1(right_degree_8[1],d,d,d);
    cipher_matrix_jiang right_degree_8_2(right_degree_8[2],d,d,d);
    cipher_matrix_jiang right_degree_8_3(right_degree_8[3],d,d,d);

    vector<Ciphertext> left_degree_16=block_256_LeftSquare(context,encoder,scale,evaluator,relin_keys,gal_keys,left_degree_8_0,left_degree_8_1,left_degree_8_2,left_degree_8_3,left_degree_8_0,left_degree_8_1,left_degree_8_2,left_degree_8_3);
    //A+...+A^16
    vector<Ciphertext> right_degree_16=block_256_LeftMulRight_AddRight(context,encoder,scale,evaluator,relin_keys,gal_keys,left_degree_8_0,left_degree_8_1,left_degree_8_2,left_degree_8_3,right_degree_8_0,right_degree_8_1,right_degree_8_2,right_degree_8_3);
    cipher_matrix_jiang left_degree_16_0(left_degree_16[0],d,d,d);
    cipher_matrix_jiang left_degree_16_1(left_degree_16[1],d,d,d);
    cipher_matrix_jiang left_degree_16_2(left_degree_16[2],d,d,d);
    cipher_matrix_jiang left_degree_16_3(left_degree_16[3],d,d,d);
    cipher_matrix_jiang right_degree_16_0(right_degree_16[0],d,d,d);
    cipher_matrix_jiang right_degree_16_1(right_degree_16[1],d,d,d);
    cipher_matrix_jiang right_degree_16_2(right_degree_16[2],d,d,d);
    cipher_matrix_jiang right_degree_16_3(right_degree_16[3],d,d,d);

    Matrix<double> computed0,computed1,computed2,computed3;
    right_degree_16_0.dec_matrix_cipher(computed0,encoder,decryptor);
    right_degree_16_1.dec_matrix_cipher(computed1,encoder,decryptor);
    right_degree_16_2.dec_matrix_cipher(computed2,encoder,decryptor);
    right_degree_16_3.dec_matrix_cipher(computed3,encoder,decryptor);
    
    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

 

    // std::ofstream outFile0("matrix_result_256_0.txt");
    // for (size_t i = 0; i < d; i++) {
    //     for (size_t j = 0; j < d; j++) {
    //         outFile0 << computed0.get(j, i) << (j < n - 1 ? " " : "");
    //     }
    //     outFile0 << std::endl;
    // }

    // outFile0.close();
    // std::ofstream outFile1("matrix_result_256_1.txt");
    // for (size_t i = 0; i < d; i++) {
    //     for (size_t j = 0; j < d; j++) {
    //         outFile1 << computed1.get(j, i) << (j < n - 1 ? " " : "");
    //     }
    //     outFile1 << std::endl;
    // }

    // outFile1.close();
    // std::ofstream outFile2("matrix_result_256_2.txt");
    // for (size_t i = 0; i < d; i++) {
    //     for (size_t j = 0; j < d; j++) {
    //         outFile2 << computed2.get(j, i) << (j < n - 1 ? " " : "");
    //     }
    //     outFile2 << std::endl;
    // }

    // outFile2.close();
    // std::ofstream outFile3("matrix_result_256_3.txt");
    // for (size_t i = 0; i < d; i++) {
    //     for (size_t j = 0; j < d; j++) {
    //         outFile3 << computed3.get(j, i) << (j < n - 1 ? " " : "");
    //     }
    //     outFile3 << std::endl;
    // }

    // outFile3.close();

    cout<<"endl"<<"---------------------------------------------------------"<<endl;
    // cout << "the true C =  A^2+A: " << endl;
    // A_true_0.print();

}
void run_pack_inversion(size_t n, size_t m, size_t p,size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    //438
//     parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30,30,30, 60}));
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 40, 40, 40, 40, 40, 40, 40, 40, 60}));
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {45, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 50}));
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 60}));
    double scale = pow(2.0, 40);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    // print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);
    cout<<"---------"<<endl;
    test_pack_inversion(context,n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}

void test_pack_inversion(SEALContext context, size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    size_t d=sqrt(encoder.slot_count());

    size_t repeat_times = 1024; // 重复次数
    std::vector<double> result_vector; // 用于存储最终的向量
    std::vector<double> result_I;
    for (size_t i = 0; i < repeat_times; ++i) {
        Matrix<double> A(n, m); // 生成矩阵 A
        Matrix<double> I(n, m);
        random_pack_matrix_generator(n, m, A); // 填充矩阵 A
        eye_matrix_generator(n,m,I);
        if(i==0){
            A.print();
        }

        // 将矩阵 A 的每一行按顺序拼接到 result_vector 中
        for (size_t row = 0; row < n; ++row) {
            for (size_t col = 0; col < m; ++col) {
                result_I.push_back(I(row,col));
                result_vector.push_back(A(row, col)); // 假设 Matrix 提供了 (row, col) 的访问方式
            }
        }
    }
    //保存明文向量
    // store_vector(result_vector);//A

    for (size_t i = 0; i < result_I.size(); ++i) {
        result_vector[i]=result_I[i]-result_vector[i];
    }
    // 输出结果向量（可选）
    cout<<"加密前: 向量前16位"<<endl;
    for (size_t i = 0; i < result_vector.size(); ++i) {
        if(i<32){
            std::cout << result_vector[i] << " ";
        }
    }
    cout<<"______________"<<endl;
    auto start_time = chrono::high_resolution_clock::now();

    seal::Plaintext plain_begin,plain_over;
    encoder.encode(result_vector,scale,plain_begin);
    seal::Ciphertext cipher_A,cipher_B;
    encryptor.encrypt(plain_begin,cipher_A);
    cipher_B=cipher_A;

    cout<<"A0"<<endl;
    // cout<<"cipher1:"<<context.get_context_data(cipher_A.parms_id())->chain_index()<<endl;
    Ciphertext cipher_A_change = pack_change_cipher_to_matrix_A(
        context,cipher_A, scale, encoder, evaluator, decryptor, gal_keys, relin_keys);
    // decrypt_and_print(cipher_A_change, decryptor, encoder);
    // cout<<"cipher1:"<<context.get_context_data(cipher_A_change.parms_id())->chain_index()<<endl;
    

    cout<<"B0"<<endl;
    // cout<<"cipher1:"<<context.get_context_data(cipher_B.parms_id())->chain_index()<<endl;
    Ciphertext cipher_B_change = pack_change_cipher_to_matrix_B(
        context,cipher_B, scale, encoder, evaluator, decryptor, gal_keys, relin_keys);
    // decrypt_and_print(cipher_B_change, decryptor, encoder);
    
    
    //A0*B0
    Ciphertext cipher_AB = pack_encrypted_matrix_mul(
        context,cipher_A_change,cipher_B_change, scale, encoder, evaluator, decryptor, gal_keys, relin_keys);
    cout<<"A^2"<<endl;
    decrypt_and_print(cipher_AB, decryptor, encoder);
    
    

    //(A+A^2)
    match_levels(context,cipher_A,cipher_AB,evaluator);
    seal:: Ciphertext cipher_A_3;
    evaluator.add(cipher_A,cipher_AB,cipher_A_3);
    cout<<"A+A^2"<<endl;
    decrypt_and_print(cipher_A_3, decryptor, encoder);

    Ciphertext cipher_AB_change = pack_change_cipher_to_matrix_A(
        context,cipher_AB, scale, encoder, evaluator, decryptor, gal_keys, relin_keys);
    Ciphertext cipher_A_3_change = pack_change_cipher_to_matrix_B(
        context,cipher_A_3, scale, encoder, evaluator, decryptor, gal_keys, relin_keys);
    //A^2(A+A^2)
    Ciphertext cipher_A_4 = pack_encrypted_matrix_mul(
        context,cipher_AB_change,cipher_A_3_change, scale, encoder, evaluator, decryptor, gal_keys, relin_keys);
    cout<<"A^3+A^4"<<endl;
    decrypt_and_print(cipher_A_4, decryptor, encoder);
    match_levels(context,cipher_A_3,cipher_A_4,evaluator);
    seal::Ciphertext cipher_result;
    evaluator.add(cipher_A_3,cipher_A_4,cipher_result);

    cout<<"cipher1:"<<context.get_context_data(cipher_result.parms_id())->chain_index()<<endl;
    decrypt_and_print(cipher_result, decryptor, encoder);


    // decryptor.decrypt(cipher_A_change,plain_over);
    // std::vector<double> result;
    // encoder.decode(plain_over,result);
    // // 输出结果向量（可选）
    // cout<<"解密后: 向量前16位"<<endl;
    // for (size_t i = 0; i < result_vector.size(); ++i) {
    //     if(i<16){
    //         std::cout << result[i] << " ";
    //     }
    //     }
    // cout<<"______________"<<endl;
    


    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "----------------- Jiang matrix multiplication --------------------" << endl;
    cout << endl
          << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "Actual Multiply:"
          << " (" << d << "x" << d << ")X(" << d << "x" << d << ")" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;
    
    cout<<d<<endl;
    std::cout << std::endl;

}


void run_dia_pack_new_encoding_inversion(size_t n, size_t m, size_t p,size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    //438
//     parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30,30,30, 60}));
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 40, 40, 40, 40, 40, 40, 40, 40, 60}));
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {45, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 50}));
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 45, 45, 45, 45, 45, 45, 45, 60}));
    double scale = pow(2.0, 45);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    // print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);
    cout<<"---------"<<endl;
    test_dia_pack_new_encoding_inversion(context,n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}


void test_dia_pack_new_encoding_inversion(SEALContext context, size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    size_t d=encoder.slot_count();
    size_t step=encoder.slot_count()/n/m;
    size_t repeat_times = d/n/m; // 重复次数
    cout<<d<<step<<n<<m<<repeat_times<<endl;
    // std::vector<double> result_vector(d); // 用于存储最终的向量
    std::vector<double> result_I(d);
    for (size_t i = 0; i < repeat_times; ++i) {
        // Matrix<double> A(n, m); // 生成矩阵 A
        Matrix<double> I(n, m);
        // random_pack_matrix_generator_new_encoding(n, m, A); // 填充矩阵 A
        eye_matrix_generator(n,m,I);
        // if(i==0){
        //     A.print();
        // }
        // 将矩阵 A 的元素按块存储到 result_vector 中
        for (size_t row = 0; row < n; ++row) {
            for (size_t col = 0; col < m; ++col) {
                size_t block_index = row * m + col; // 当前块的索引
                result_I[block_index * repeat_times + i] = I.get(row, col);
                // result_vector[block_index * repeat_times + i] = A.get(row, col);
            }
        }
    }
    //保存明文向量
    // store_vector_new_encoding(result_vector);//A
    // 计算I-A
    // for (size_t i = 0; i < result_I.size(); ++i) {
    //     result_vector[i]=result_I[i]-result_vector[i];
    // }
    // // 输出结果向量（可选）
    // cout<<"加密前: 每块第一位:"<<endl;
    // for (size_t i = 0; i < result_vector.size(); ++i) {
    //     if(i%step==0)
    //     {
    //         std::cout << result_vector[i] << " ";
    //     }
    // }
    std::vector<double> result_vector; 
    try {
        // 文件的相对路径
        // std::string file_path = "../build/dia_input/dia_input_c_(128,128).txt";
        std::string file_path = "../build/dia_triangle/input_dia_(4,4).txt";

        // 调用函数读取数据
        result_vector = read_data_from_file(file_path);

        std::cout << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "发生错误: " << e.what() << std::endl;
    }


    cout<<"______________"<<endl;
    auto start_time = chrono::high_resolution_clock::now();

    seal::Plaintext plain_begin,plain_over;
    encoder.encode(result_vector,scale,plain_begin);
    seal::Ciphertext cipher_A,cipher_B;
    encryptor.encrypt(plain_begin,cipher_A);
    cipher_B=cipher_A;
    cout<<"cipher1:"<<context.get_context_data(cipher_A.parms_id())->chain_index()<<endl;
    // cout<<"A0"<<endl;
    Ciphertext cipher_A_change = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_A, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_A_change, decryptor, encoder,n);

    // cout<<"B0"<<endl;
    Ciphertext cipher_B_change = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_B, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_B_change, decryptor, encoder,n);
    
    //A0*B0
    Ciphertext cipher_AB = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_A_change,cipher_B_change, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // cout<<"A^2"<<endl;
    // decrypt_and_print_new_encoding(cipher_AB, decryptor, encoder,n);
    
    //(A+A^2)
    match_levels(context,cipher_A,cipher_AB,evaluator);
    seal:: Ciphertext cipher_A_3;
    evaluator.add(cipher_A,cipher_AB,cipher_A_3);
    // cout<<"A+A^2"<<endl;
    // decrypt_and_print_new_encoding(cipher_A_3, decryptor, encoder);

   Ciphertext cipher_AB_change = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_AB, scale, encoder, evaluator, decryptor, gal_keys,  relin_keys,n);
    Ciphertext cipher_A_3_change = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_A_3, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    //A^2(A+A^2)
    Ciphertext cipher_A_4 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_AB_change,cipher_A_3_change, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // cout<<"A^3+A^4"<<endl;
    // decrypt_and_print_new_encoding(cipher_A_4, decryptor, encoder);
    match_levels(context,cipher_A_3,cipher_A_4,evaluator);
    seal::Ciphertext cipher_result;
    evaluator.add(cipher_A_3,cipher_A_4,cipher_result);

    cout<<"cipher1:"<<context.get_context_data(cipher_result.parms_id())->chain_index()<<endl;
    decrypt_and_print_new_encoding_over(cipher_result, decryptor, encoder,result_I,n);







    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "----------------- Jiang matrix multiplication --------------------" << endl;
    cout << endl
          << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "Actual Multiply:"
          << " (" << d << "x" << d << ")X(" << d << "x" << d << ")" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;
    
    std::cout << std::endl;

}
void run_real_pack_new_encoding_inversion(size_t n, size_t m, size_t p,size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    //438
//     parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30,30,30, 60}));
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 40, 40, 40, 40, 40, 40, 40, 40, 60}));
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {45, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 50}));
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 60}));
    double scale = pow(2.0, 45);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    // print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);
    cout<<"---------"<<endl;
    test_real_pack_new_encoding_inversion(context,n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}


void test_real_pack_new_encoding_inversion(SEALContext context, size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    size_t d=encoder.slot_count();
    size_t step=encoder.slot_count()/n/m;
    size_t repeat_times = d/n/m; // 重复次数
    cout<<d<<step<<n<<m<<repeat_times<<endl;
    // std::vector<double> result_vector(d); // 用于存储最终的向量
    std::vector<double> result_I(d);
    for (size_t i = 0; i < repeat_times; ++i) {
        Matrix<double> A(n, m); // 生成矩阵 A
        Matrix<double> I(n, m);
        // random_pack_matrix_generator_new_encoding(n, m, A); // 填充矩阵 A
        eye_matrix_generator(n,m,I);
        // if(i==0){
        //     A.print();
        // }
        // 将矩阵 A 的元素按块存储到 result_vector 中
        for (size_t row = 0; row < n; ++row) {
            for (size_t col = 0; col < m; ++col) {
                size_t block_index = row * m + col; // 当前块的索引
                result_I[block_index * repeat_times + i] = I.get(row, col);
                // result_vector[block_index * repeat_times + i] = A.get(row, col);
            }
        }
    }
    //保存明文向量
    // store_vector_new_encoding(result_vector);//A
    //计算I-A
    // for (size_t i = 0; i < result_I.size(); ++i) {
    //     result_vector[i]=result_I[i]-result_vector[i];
    // }
    // // 输出结果向量（可选）
    // cout<<"加密前: 每块第一位:"<<endl;
    // for (size_t i = 0; i < result_vector.size(); ++i) {
    //     if(i%step==0)
    //     {
    //         std::cout << result_vector[i] << " ";
    //     }
    // }
    std::vector<double> result_vector; 
    try {
        // 文件的相对路径
        // std::string file_path = "../build/real_input/real_input_c_(4,4).txt";//random
        std::string file_path = "../build/dia_movie/dia_movie_c_all_852_(4,4).txt";

        // 调用函数读取数据
        result_vector = read_data_from_file(file_path);

        std::cout << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "发生错误: " << e.what() << std::endl;
    }


    // cout<<"______________"<<endl;
    auto start_time = chrono::high_resolution_clock::now();

    seal::Plaintext plain_begin,plain_over;
    encoder.encode(result_vector,scale,plain_begin);
    seal::Ciphertext cipher_A,cipher_B;
    encryptor.encrypt(plain_begin,cipher_A);
    cipher_B=cipher_A;
    cout<<"cipher1:"<<context.get_context_data(cipher_A.parms_id())->chain_index()<<endl;
    // cout<<"A0"<<endl;
    Ciphertext cipher_A_change = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_A, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_A_change, decryptor, encoder,n);

    // cout<<"B0"<<endl;
    Ciphertext cipher_B_change = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_B, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_B_change, decryptor, encoder,n);
    
    //A0*B0
    Ciphertext cipher_AB = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_A_change,cipher_B_change, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // cout<<"A^2"<<endl;
    // decrypt_and_print_new_encoding(cipher_AB, decryptor, encoder,n)
    //A^4
    Ciphertext cipher_AB_change_1=pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_AB, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    Ciphertext cipher_AB_change_2=pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_AB, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    Ciphertext cipher_A_power_4 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_AB_change_1,cipher_AB_change_2, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    //(A+A^2)
    match_levels(context,cipher_A,cipher_AB,evaluator);
    seal:: Ciphertext cipher_A_3;
    evaluator.add(cipher_A,cipher_AB,cipher_A_3);
    // cout<<"A+A^2"<<endl;
    // decrypt_and_print_new_encoding(cipher_A_3, decryptor, encoder);

    Ciphertext cipher_AB_change = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_AB, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    Ciphertext cipher_A_3_change = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_A_3, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    //A^2(A+A^2)

    Ciphertext cipher_A_4 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_AB_change,cipher_A_3_change, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // cout<<"A^3+A^4"<<endl;
    // decrypt_and_print_new_encoding(cipher_A_4, decryptor, encoder);
    match_levels(context,cipher_A_3,cipher_A_4,evaluator);
    seal::Ciphertext cipher_result;
    evaluator.add(cipher_A_3,cipher_A_4,cipher_result);

    // cout<<"cipher1:"<<context.get_context_data(cipher_result.parms_id())->chain_index()<<endl;
    // decrypt_and_print_new_encoding_over(cipher_result, decryptor, encoder,result_I,n);


    seal:: Plaintext plain_A1,plain_Y1,plain_A1_1,plain_Y1_1;
    seal:: Ciphertext cipher_A1_tmp,cipher_Y1_tmp;
    decryptor.decrypt(cipher_A_power_4,plain_A1);
    decryptor.decrypt(cipher_result,plain_Y1);
    std::vector<double> result_A1;
    std::vector<double> result_Y1;
    encoder.decode(plain_A1, result_A1);
    encoder.decode(plain_Y1, result_Y1);

   
   

    seal:: Plaintext plain_left_4,plain_right_4;
    seal:: Ciphertext cipher_left_4,cipher_right_4;
    cipher_left_4=cipher_A_power_4;
    cipher_right_4=cipher_result;


    Ciphertext cipher_left_4_change1 = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_left_4, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    Ciphertext cipher_left_4_change2 = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_left_4, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    Ciphertext cipher_right_4_change = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_right_4, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    //A^5+..+A^8
    Ciphertext cipher_right_8 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_left_4_change1,cipher_right_4_change, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    //A^8
    Ciphertext cipher_left_8 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_left_4_change1,cipher_left_4_change2, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    match_levels(context,cipher_right_4,cipher_right_8,evaluator);
    //A^1+..+A^8
    seal:: Ciphertext cipher_A_1_8;
    evaluator.add(cipher_right_4,cipher_right_8,cipher_A_1_8);

    // decrypt_and_print_new_encoding_over(cipher_A_1_8, decryptor, encoder,result_I,n);

    // auto end_time = chrono::high_resolution_clock::now();
    // auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);

  
    // cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;
    
    // std::cout << std::endl;

    Ciphertext cipher_left_8_change =pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_left_8, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    Ciphertext cipher_A_1_8_change = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_A_1_8, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    //A^9+...+A^16
    Ciphertext cipher_right_16 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_left_8_change,cipher_A_1_8_change, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    match_levels(context,cipher_A_1_8,cipher_right_16,evaluator);
    seal::Ciphertext fin_result;
    evaluator.add(cipher_A_1_8,cipher_right_16,fin_result);

    cout<<"cipher1:"<<context.get_context_data(fin_result.parms_id())->chain_index()<<endl;
    
    decrypt_and_print_new_encoding_over(fin_result, decryptor, encoder,result_I,n);

    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);

  
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;
    
    std::cout << std::endl;

}
void run_triangle_inversion(size_t n, size_t m, size_t p,size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    //438
//     parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30,30,30, 60}));
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 40, 40, 40, 40, 40, 40, 40, 40, 60}));
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {45, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 50}));
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 45, 45, 45, 45, 45, 45, 45, 60}));
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 45, 45, 45, 45, 45, 60}));
    double scale = pow(2.0, 45);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    // print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);
    cout<<"---------"<<endl;
    test_triangle_inversion(context,n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}
void test_triangle_inversion(SEALContext context, size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    size_t d=encoder.slot_count();
    cout<<d<<endl;
    // size_t repeat_times=d/n/m;
    // std::vector<double> result_vector; // 用于存储最终的向量
    // for (size_t i = 0; i < repeat_times; ++i) {
    //     Matrix<double> A(n, m); // 生成矩阵 A
    //     random_pack_triangle_matrix_generator(n, m, A); // 填充矩阵 A
    //     if(i==0){
    //         A.print();
    //     }

    //     // 将矩阵 A 的每一行按顺序拼接到 result_vector 中
    //     for (size_t row = 0; row < n; ++row) {
    //         for (size_t col = 0; col < m; ++col) {
    //             result_vector.push_back(A(row, col)); // 假设 Matrix 提供了 (row, col) 的访问方式
    //         }
    //     }
    // }
    std::vector<double> result_vector; 
    try {
        // 文件的相对路径
        // std::string file_path = "../build/dia_input/dia_input_c_(128,128).txt";
        std::string file_path = "../build/dia_triangle/input_triangle_(4,4).txt";

        // 调用函数读取数据
        result_vector = read_data_from_file(file_path);

        std::cout << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "发生错误: " << e.what() << std::endl;
    }
    //保存明文向量
    // store_vector_triangle(result_vector);//A
    auto start_time = chrono::high_resolution_clock::now();
    seal::Plaintext plain_begin,p1,p2,p3,p4;
    encoder.encode(result_vector,scale,plain_begin);
    seal::Ciphertext cipher_A,cipher_1,cipher_2,cipher_3,cipher_4;
    encryptor.encrypt(plain_begin,cipher_A);
    // cout<<"cipherA:"<<context.get_context_data(cipher_A.parms_id())->chain_index()<<endl;
    // decrypt_and_print(cipher_A, decryptor, encoder);

    vector<double> e1(d);
    for(size_t i =0;i<d;i++){
        if(i%16==0){
            e1[i]=1;
        }
        else{
            e1[i]=0;
        }
    }
    vector<double> e2(d);
    for(size_t i =0;i<d;i++){
        if(i%16==5){
            e2[i]=1;
        }
        else{
            e2[i]=0;
        }
    }
    vector<double> e3(d);
    for(size_t i =0;i<d;i++){
        if(i%16==10){
            e3[i]=1;
        }
        else{
            e3[i]=0;
        }
    }
    vector<double> e4(d);
    for(size_t i =0;i<d;i++){
        if(i%16==15){
            e4[i]=1;
        }
        else{
            e4[i]=0;
        }
    }
    // encoder.encode(e2,scale,p2);
    // encoder.encode(e3,scale,p3);
    // encoder.encode(e4,scale,p4);
    //p1

    

    //p2
    seal::Ciphertext cipher_A_21;
    seal::Plaintext p21;
    vector<double> e21(d);
    for(size_t i =0;i<d;i++){
        if(i%16==4){
            e21[i]=1;
        }
        else{
            e21[i]=0;
        }
    }
    // for(size_t i =0 ;i<32;i++)
    // {
    //     cout<<e21[i]<<" ";
    //     if(i%16==15){
    //         cout<<endl;
    //     }
    // }
    encoder.encode(e21,cipher_A.scale(),p21);
    evaluator.mod_switch_to_inplace(p21,cipher_A.parms_id());
    evaluator.multiply_plain(cipher_A,p21,cipher_A_21);
    evaluator.rescale_to_next_inplace(cipher_A_21);
    encoder.encode(e2,cipher_A_21.scale(),p2);
    evaluator.mod_switch_to_inplace(p2,cipher_A_21.parms_id());
 
    evaluator.add_plain(cipher_A_21,p2,cipher_2);
    // decrypt_and_print(cipher_2, decryptor, encoder);
    // cout<<"cipher2:"<<context.get_context_data(cipher_2.parms_id())->chain_index()<<endl;
    //p3
    seal::Ciphertext cipher_A_31,cipher_A_32;
    seal::Plaintext p31,p32;
    vector<double> e31(d),e32(d);
    for(size_t i =0;i<d;i++){
        if(i%16==8){
            e31[i]=1;
        }
        else{
            e31[i]=0;
        }
    }
    for(size_t i =0;i<d;i++){
        if(i%16==9){
            e32[i]=1;
        }
        else{
            e32[i]=0;
        }
    }
    encoder.encode(e31,cipher_A.scale(),p31);
    evaluator.mod_switch_to_inplace(p31,cipher_A.parms_id());
    evaluator.multiply_plain(cipher_A,p31,cipher_A_31);
    evaluator.rescale_to_next_inplace(cipher_A_31);
    // decrypt_and_print(cipher_A_31, decryptor, encoder);

    encoder.encode(e32,cipher_A.scale(),p32);
    evaluator.mod_switch_to_inplace(p32,cipher_A.parms_id());
    evaluator.multiply_plain(cipher_A,p32,cipher_A_32);
    evaluator.rescale_to_next_inplace(cipher_A_32);
    // decrypt_and_print(cipher_A_32, decryptor, encoder);
    seal::Ciphertext tmp31,tmp32;
    tmp31=cipher_A_31;
    tmp32=cipher_A_32;

    cipher_3=rotate_as_16(cipher_A_21, -n, scale, encoder, evaluator, decryptor,gal_keys, relin_keys);
    seal:: Ciphertext tmp;
    tmp = rotate_as_16(cipher_A_32, 1, scale, encoder, evaluator, decryptor,gal_keys, relin_keys);
    evaluator.multiply_inplace(cipher_3,tmp);
    evaluator.relinearize_inplace(cipher_3, relin_keys);
    evaluator.rescale_to_next_inplace(cipher_3);
    
    encoder.encode(e3,cipher_3.scale(),p3);
    evaluator.mod_switch_to_inplace(p3,cipher_3.parms_id());
    evaluator.add_plain_inplace(cipher_3,p3);
    // decrypt_and_print(cipher_3, decryptor, encoder);

    match_levels(context,cipher_A_32,cipher_3,evaluator);
    evaluator.add_inplace(cipher_3,cipher_A_32);

    match_levels(context,cipher_A_31,cipher_3,evaluator);
    evaluator.add_inplace(cipher_3,cipher_A_31);

    // cout<<"cipher3:"<<context.get_context_data(cipher_3.parms_id())->chain_index()<<endl;
    // decrypt_and_print(cipher_3, decryptor, encoder);

    //p4
    seal::Ciphertext cipher_A_41,cipher_A_42,cipher_A_43;
    seal::Plaintext p41,p42,p43;
    vector<double> e41(d),e42(d),e43(d);
    for(size_t i =0;i<d;i++){
        if(i%16==12){
            e41[i]=1;
        }
        else{
            e41[i]=0;
        }
    }
    for(size_t i =0;i<d;i++){
        if(i%16==13){
            e42[i]=1;
        }
        else{
            e42[i]=0;
        }
    }
    for(size_t i =0;i<d;i++){
        if(i%16==14){
            e43[i]=1;
        }
        else{
            e43[i]=0;
        }
    }
    
    encoder.encode(e41,cipher_A.scale(),p41);
    evaluator.mod_switch_to_inplace(p41,cipher_A.parms_id());
    evaluator.multiply_plain(cipher_A,p41,cipher_A_41);
    evaluator.rescale_to_next_inplace(cipher_A_41);
    // decrypt_and_print(cipher_A_41, decryptor, encoder);

    encoder.encode(e42,cipher_A.scale(),p42);
    evaluator.mod_switch_to_inplace(p42,cipher_A.parms_id());
    evaluator.multiply_plain(cipher_A,p42,cipher_A_42);
    evaluator.rescale_to_next_inplace(cipher_A_42);
    // decrypt_and_print(cipher_A_42,decryptor,encoder);

    encoder.encode(e43,cipher_A.scale(),p43);
    evaluator.mod_switch_to_inplace(p43,cipher_A.parms_id());
    evaluator.multiply_plain(cipher_A,p43,cipher_A_43);
    evaluator.rescale_to_next_inplace(cipher_A_43);
    // decrypt_and_print(cipher_A_43,decryptor,encoder);
    seal::Ciphertext tmp42,tmp43;
    tmp42=cipher_A_42;
    tmp43=cipher_A_43;

    seal::Ciphertext tmp_43,tmp_32,tmp_21,tmp_42,tmp_31;
    tmp_21=rotate_as_16(cipher_A_21, -2*n, scale, encoder, evaluator, decryptor,gal_keys, relin_keys);
    tmp_32=rotate_as_16(tmp32, -3, scale, encoder, evaluator, decryptor,gal_keys, relin_keys);
    tmp_43=rotate_as_16(cipher_A_43, 2, scale, encoder, evaluator, decryptor,gal_keys, relin_keys);
    evaluator.multiply(tmp_43,tmp_32,cipher_4);
    evaluator.relinearize_inplace(cipher_4, relin_keys);
    evaluator.rescale_to_next_inplace(cipher_4);
    match_levels(context,tmp_21,cipher_4,evaluator);
    evaluator.multiply_inplace(cipher_4,tmp_21);
    evaluator.relinearize_inplace(cipher_4, relin_keys);
    evaluator.rescale_to_next_inplace(cipher_4);
    //e4
    encoder.encode(e4,cipher_4.scale(),p4);
    evaluator.mod_switch_to_inplace(p4,cipher_4.parms_id());
    evaluator.add_plain_inplace(cipher_4,p4);
    //A41
    match_levels(context,cipher_A_41,cipher_4,evaluator);
    evaluator.add_inplace(cipher_4,cipher_A_41);
    //A42 e2
    match_levels(context,cipher_A_42,cipher_4,evaluator);
    evaluator.add_inplace(cipher_4,cipher_A_42);
    //A42 A21 e1
    tmp_21 = rotate_as_16(cipher_A_21, -8, scale, encoder, evaluator, decryptor,gal_keys, relin_keys);
    tmp_42 = rotate_as_16(tmp42, 1, scale, encoder, evaluator, decryptor,gal_keys, relin_keys);
    evaluator.multiply_inplace(tmp_42,tmp_21);
    evaluator.relinearize_inplace(tmp_42, relin_keys);
    evaluator.rescale_to_next_inplace(tmp_42);
    match_levels(context,tmp_42,cipher_4,evaluator);
    evaluator.add_inplace(cipher_4,tmp_42);

    //A43 e3
    match_levels(context,cipher_A_43,cipher_4,evaluator);
    evaluator.add_inplace(cipher_4,cipher_A_43);

    //A43 a31 e1
    tmp_31=rotate_as_16(tmp31, -n, scale, encoder, evaluator, decryptor,gal_keys, relin_keys);
    evaluator.multiply_inplace(tmp_31,tmp_43);
    evaluator.relinearize_inplace(tmp_31, relin_keys);
    evaluator.rescale_to_next_inplace(tmp_31);
    match_levels(context,tmp_31,cipher_4,evaluator);
    evaluator.add_inplace(cipher_4,tmp_31);
    //A43 A32 e2
    tmp_43=rotate_as_16(tmp43, 1, scale, encoder, evaluator, decryptor,gal_keys, relin_keys);
    tmp_32=rotate_as_16(tmp32, -4, scale, encoder, evaluator, decryptor,gal_keys, relin_keys);
    evaluator.multiply_inplace(tmp_43,tmp_32);
    evaluator.relinearize_inplace(tmp_43, relin_keys);
    evaluator.rescale_to_next_inplace(tmp_43);
    match_levels(context,tmp_43,cipher_4,evaluator);
    evaluator.add_inplace(cipher_4,tmp_43);

    //all add
    encoder.encode(e1,cipher_4.scale(),p1);
    evaluator.mod_switch_to_inplace(p1,cipher_4.parms_id());
    evaluator.add_plain_inplace(cipher_4,p1);

    match_levels(context,cipher_2,cipher_4,evaluator);
    evaluator.add_inplace(cipher_4,cipher_2);

    match_levels(context,cipher_3,cipher_4,evaluator);
    evaluator.add_inplace(cipher_4,cipher_3);


    // cout<<"cipher4:"<<context.get_context_data(cipher_4.parms_id())->chain_index()<<endl;












    decrypt_and_print(cipher_4, decryptor, encoder);
    // decrypt_and_print(cipher_A, decryptor, encoder);

    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << endl
         << "----------------- Jiang matrix multiplication --------------------" << endl;
    cout << endl
          << "(n, m, p): (" << n << ", " << m << ", " << p << ")" << endl;
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;
    
    std::cout << std::endl;

}

void run_tae_inversion(size_t n, size_t m, size_t p,size_t poly_modulus_degree)
{

    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    //438
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 60}));
    double scale = pow(2.0, 45);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    // print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);
    cout<<"---------"<<endl;
    test_tae_inversion(context,n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}


void test_tae_inversion(SEALContext context, size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    size_t d=encoder.slot_count();
    size_t step=encoder.slot_count()/n/m;
    size_t repeat_times = d/n/m; // 重复次数
    // cout<<d<<n<<m<<repeat_times<<endl;
    // std::vector<double> result_vector(d); // 用于存储最终的向量
    std::vector<double> result_I(d);
    for (size_t i = 0; i < repeat_times; ++i) {
        Matrix<double> A(n, m); // 生成矩阵 A
        Matrix<double> I(n, m);
        // random_pack_matrix_generator_new_encoding(n, m, A); // 填充矩阵 A
        // random_tae_matrix_generator(n,m,A);
        eye_matrix_generator(n,m,I);
        // if(i==0){
        //     A.print();
        // }
        

        // 将矩阵 A 的元素按块存储到 result_vector 中
        for (size_t row = 0; row < n; ++row) {
            for (size_t col = 0; col < m; ++col) {
                size_t block_index = row * m + col; // 当前块的索引
                result_I[block_index * repeat_times + i] = I.get(row, col);
                // result_vector[block_index * repeat_times + i] = A.get(row, col);
            }
        }
    }
    //保存明文向量
    // store_vector_new_encoding(result_vector);//A
    std::vector<double> result_vector; // 用于存储最终的向量
    try {
        // 文件的相对路径
        std::string file_path = "../build/dia_triangle/input_tae_(4,4).txt";

        // 调用函数读取数据
        result_vector = read_data_from_file(file_path);

        std::cout << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "发生错误: " << e.what() << std::endl;
    }

    std::vector<double> Y0_data;
    try {
        // 文件的相对路径
        std::string file_path = "../build/dia_triangle/input_tae_Y0_(4,4).txt";

        // 调用函数读取数据
        Y0_data = read_data_from_file(file_path);

        std::cout << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "发生错误: " << e.what() << std::endl;
    }


    // 输出结果向量（可选）
    cout<<"加密前: 每块第一位:"<<endl;
    for (size_t i = 0; i < result_vector.size(); ++i) {
        if(i%step==0)
        {
            std::cout << result_vector[i] << " ";
        }
    }
    cout<<Y0_data.size()<<endl;
    
    // cout<<"______________"<<endl;
    auto start_time0 = chrono::high_resolution_clock::now();

    seal::Plaintext plain_Y0,plain_I,plain_A,plain_0;
    encoder.encode(result_vector,scale,plain_A);
    encoder.encode(Y0_data,scale,plain_Y0);
    
    
    seal::Ciphertext cipher_A,cipher_Y0;
    encryptor.encrypt(plain_Y0,cipher_Y0);
    encryptor.encrypt(plain_A,cipher_A);
    cout<<"cipher1:"<<context.get_context_data(cipher_A.parms_id())->chain_index()<<endl;
    // decrypt_and_print_new_encoding(cipher_Y0, decryptor, encoder,n);


    Ciphertext cipher_A_change = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_A, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_A_change, decryptor, encoder,n);


    Ciphertext cipher_Y0_tmp = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_Y0, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_B_change, decryptor, encoder,n);
    
    //A0=A*Y0
    Ciphertext cipher_A0 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_A_change,cipher_Y0_tmp, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_A0, decryptor, encoder,n);

    encoder.encode(-1,cipher_A0.scale(),plain_0);
    evaluator.mod_switch_to_inplace(plain_0,cipher_A0.parms_id());
    evaluator.multiply_plain_inplace(cipher_A0,plain_0);


    encoder.encode(result_I,cipher_A0.scale(),plain_I);
    evaluator.mod_switch_to_inplace(plain_I,cipher_A0.parms_id());
    evaluator.add_plain_inplace(cipher_A0,plain_I);
    evaluator.rescale_to_next_inplace(cipher_A0);
    // decrypt_and_print_new_encoding(cipher_A0, decryptor, encoder,n);

    

    match_levels(context,cipher_Y0,cipher_A0,evaluator);
    // decrypt_and_print_new_encoding(cipher_Y0, decryptor, encoder,n);

    Ciphertext cipher_Y0_change = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_Y0, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_A_change, decryptor, encoder,n);

    seal::Plaintext plain_I_0;
    seal::Ciphertext cipher_A0_I;
    encoder.encode(result_I,cipher_A0.scale(),plain_I_0);
    evaluator.mod_switch_to_inplace(plain_I_0,cipher_A0.parms_id());
    evaluator.add_plain(cipher_A0,plain_I_0,cipher_A0_I);
    // decrypt_and_print_new_encoding(cipher_A0_I, decryptor, encoder,n);

    Ciphertext cipher_A0_change_I = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_A0_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_B_change, decryptor, encoder,n);
    
    //Y1=Y0*(I+A0)
    Ciphertext cipher_Y1 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_Y0_change,cipher_A0_change_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_Y1, decryptor, encoder,n);


    //A1=A0^2
    seal::Ciphertext cipher_B;
    cipher_B=cipher_A0;
    Ciphertext cipher_B_change = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_B, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    Ciphertext cipher_A0_change = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_A0, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    Ciphertext cipher_A1 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_B_change,cipher_A0_change, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
   

    // cout<<"cipher1:"<<context.get_context_data(cipher_A1.parms_id())->chain_index()<<endl;
    // cout<<"cipher1:"<<context.get_context_data(cipher_Y1.parms_id())->chain_index()<<endl;
    
    // decrypt_and_print_new_encoding(cipher_A1, decryptor, encoder,n);
    // decrypt_and_print_new_encoding(cipher_Y1, decryptor, encoder,n);


  




    // cout<<"cipher1:"<<context.get_context_data(cipher_A1_tmp.parms_id())->chain_index()<<endl;
    // decrypt_and_print_new_encoding(cipher_A1, decryptor, encoder,n);


    seal::Plaintext plain_I_1;
    seal::Ciphertext cipher_A1_I;
    encoder.encode(result_I,cipher_A1.scale(),plain_I_1);
    evaluator.mod_switch_to_inplace(plain_I_1,cipher_A1.parms_id());
    evaluator.add_plain(cipher_A1,plain_I_1,cipher_A1_I);
    // decrypt_and_print_new_encoding(cipher_A0_I, decryptor, encoder,n);
    Ciphertext cipher_Y1_tmp_change = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_Y1, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_A_change, decryptor, encoder,n);


    Ciphertext cipher_A1_change_I = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_A1_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_B_change, decryptor, encoder,n);
    
    //Y1=Y0*(I+A0)
    Ciphertext cipher_Y2 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_Y1_tmp_change,cipher_A1_change_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);

    //A1=A0^2

    Ciphertext cipher_B_change_1 = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_A1, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    Ciphertext cipher_A1_change_1 = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_A1, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    Ciphertext cipher_A2 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_B_change_1,cipher_A1_change_1, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
   
    // cout<<"cipher1:"<<context.get_context_data(cipher_A2.parms_id())->chain_index()<<endl;
    // cout<<"cipher1:"<<context.get_context_data(cipher_Y2.parms_id())->chain_index()<<endl;
    
    // decrypt_and_print_new_encoding(cipher_A2, decryptor, encoder,n);
    // decrypt_and_print_new_encoding(cipher_Y2, decryptor, encoder,n);

    seal::Plaintext plain_I_2;
    seal::Ciphertext cipher_A2_I;
    encoder.encode(result_I,cipher_A2.scale(),plain_I_2);
    evaluator.mod_switch_to_inplace(plain_I_2,cipher_A2.parms_id());
    evaluator.add_plain(cipher_A2,plain_I_2,cipher_A2_I);

    Ciphertext cipher_Y2_change = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_Y2, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_A_change, decryptor, encoder,n);


    Ciphertext cipher_A2_change_I = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_A2_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_B_change, decryptor, encoder,n);
    
    //Y1=Y0*(I+A0)
    Ciphertext cipher_Y3 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_Y2_change,cipher_A2_change_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);

    Ciphertext cipher_A2_change_1 = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_A2, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    Ciphertext cipher_A2_change_2 = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_A2, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    Ciphertext cipher_A3 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_A2_change_1,cipher_A2_change_2, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
   
    // cout<<"cipher1:"<<context.get_context_data(cipher_A3.parms_id())->chain_index()<<endl;
    // cout<<"cipher1:"<<context.get_context_data(cipher_Y3.parms_id())->chain_index()<<endl;
    
    // decrypt_and_print_new_encoding(cipher_A3, decryptor, encoder,n);
    // decrypt_and_print_new_encoding(cipher_Y3, decryptor, encoder,n);





    seal::Plaintext plain_I_3;
    seal::Ciphertext cipher_A3_I;
    encoder.encode(result_I,cipher_A3.scale(),plain_I_3);
    evaluator.mod_switch_to_inplace(plain_I_3,cipher_A3.parms_id());
    evaluator.add_plain(cipher_A3,plain_I_3,cipher_A3_I);
    // decrypt_and_print_new_encoding(cipher_A0_I, decryptor, encoder,n);
    Ciphertext cipher_Y3_tmp_change = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_Y3, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_A_change, decryptor, encoder,n);


    Ciphertext cipher_A3_change_I = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_A3_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_B_change, decryptor, encoder,n);
    
    //Y1=Y0*(I+A0)
    Ciphertext cipher_Y4 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_Y3_tmp_change,cipher_A3_change_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);

    //A1=A0^2

    Ciphertext cipher_A3_change_1 = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_A3, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    Ciphertext cipher_A3_change_2 = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_A3, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    Ciphertext cipher_A4 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_A3_change_1,cipher_A3_change_2, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    

    auto end_time0 = chrono::high_resolution_clock::now();
    auto time_diff0 = chrono::duration_cast<chrono::microseconds>(end_time0 - start_time0);


    // decrypt_and_print_new_encoding(cipher_A4, decryptor, encoder,n);
    // decrypt_and_print_new_encoding(cipher_Y4, decryptor, encoder,n);
    seal:: Plaintext plain_A4,plain_Y4,plain_A4_1,plain_Y4_1;
    seal:: Ciphertext cipher_A4_tmp,cipher_Y4_tmp;
    decryptor.decrypt(cipher_A4,plain_A4);
    decryptor.decrypt(cipher_Y4,plain_Y4);

    std::vector<double> result_A4;
    std::vector<double> result_Y4;
    encoder.decode(plain_A4, result_A4);
    encoder.decode(plain_Y4, result_Y4);

    auto start_time1 = chrono::high_resolution_clock::now();

    encoder.encode(result_A4,scale,plain_A4_1);
    encoder.encode(result_Y4,scale,plain_Y4_1);

    encryptor.encrypt(plain_A4_1,cipher_A4_tmp);
    encryptor.encrypt(plain_Y4_1,cipher_Y4_tmp);


    seal::Plaintext plain_I_4;
    seal::Ciphertext cipher_A4_I;
    encoder.encode(result_I,cipher_A4_tmp.scale(),plain_I_4);
    evaluator.mod_switch_to_inplace(plain_I_4,cipher_A4_tmp.parms_id());
    evaluator.add_plain(cipher_A4_tmp,plain_I_4,cipher_A4_I);

    Ciphertext cipher_Y4_change = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_Y4_tmp, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_A_change, decryptor, encoder,n);


    Ciphertext cipher_A4_change_I = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_A4_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_B_change, decryptor, encoder,n);
    
    //Y1=Y0*(I+A0)
    Ciphertext cipher_Y5 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_Y4_change,cipher_A4_change_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);

    Ciphertext cipher_A4_change_1 = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_A4_tmp, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    Ciphertext cipher_A4_change_2 = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_A4_tmp, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    Ciphertext cipher_A5 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_A4_change_1,cipher_A4_change_2, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
   
    // cout<<"cipher1:"<<context.get_context_data(cipher_A5.parms_id())->chain_index()<<endl;
    // cout<<"cipher1:"<<context.get_context_data(cipher_Y5.parms_id())->chain_index()<<endl;
    
    // decrypt_and_print_new_encoding(cipher_A5, decryptor, encoder,n);
    // decrypt_and_print_new_encoding(cipher_Y5, decryptor, encoder,n);




    // //到这里是r=5


    // auto time_all = time_diff0+time_diff1+time_diff2;
    // cout <<"the all time is :"<<time_all.count() / 1e6 << " s" << endl;
    



    seal::Plaintext plain_I_5;
    seal::Ciphertext cipher_A5_I;
    encoder.encode(result_I,cipher_A5.scale(),plain_I_5);
    evaluator.mod_switch_to_inplace(plain_I_5,cipher_A5.parms_id());
    evaluator.add_plain(cipher_A5,plain_I_5,cipher_A5_I);
    // decrypt_and_print_new_encoding(cipher_A0_I, decryptor, encoder,n);
    Ciphertext cipher_Y5_tmp_change = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_Y5, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_A_change, decryptor, encoder,n);
    Ciphertext cipher_A5_change_I = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_A5_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_B_change, decryptor, encoder,n);
    
    //Y1=Y0*(I+A0)
    Ciphertext cipher_Y6 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_Y5_tmp_change,cipher_A5_change_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);

    Ciphertext cipher_A5_change_1 = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_A5, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    Ciphertext cipher_A5_change_2 = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_A5, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    Ciphertext cipher_A6 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_A5_change_1,cipher_A5_change_2, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    seal::Plaintext plain_I_6;
    seal::Ciphertext cipher_A6_I;
    encoder.encode(result_I,cipher_A6.scale(),plain_I_6);
    evaluator.mod_switch_to_inplace(plain_I_6,cipher_A6.parms_id());
    evaluator.add_plain(cipher_A6,plain_I_6,cipher_A6_I);

    Ciphertext cipher_Y6_change = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_Y6, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_A_change, decryptor, encoder,n);


    Ciphertext cipher_A6_change_I = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_A6_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // decrypt_and_print_new_encoding(cipher_B_change, decryptor, encoder,n);
    Ciphertext cipher_Y7 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_Y6_change,cipher_A6_change_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);

    Ciphertext cipher_A6_change_1 = pack_change_cipher_to_matrix_A_new_encoding(
        context,cipher_A6, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    Ciphertext cipher_A6_change_2 = pack_change_cipher_to_matrix_B_new_encoding(
        context,cipher_A6, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    Ciphertext cipher_A7 = pack_encrypted_matrix_mul_new_encoding(
        context,cipher_A6_change_1,cipher_A6_change_2, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    

    decrypt_and_print_new_encoding(cipher_Y7, decryptor, encoder,n);

    //到这里是r=7
    auto end_time1 = chrono::high_resolution_clock::now();
    auto time_diff1 = chrono::duration_cast<chrono::microseconds>(end_time1 - start_time1);


    auto time_all = time_diff0+time_diff1;
    cout <<"the all time is :"<<time_all.count() / 1e6 << " s" << endl;
    std::cout << std::endl;


    // seal::Plaintext plain_I_7;
    // seal::Ciphertext cipher_A7_I;
    // encoder.encode(result_I,cipher_A7.scale(),plain_I_7);
    // evaluator.mod_switch_to_inplace(plain_I_7,cipher_A7.parms_id());
    // evaluator.add_plain(cipher_A7,plain_I_7,cipher_A7_I);
    // // decrypt_and_print_new_encoding(cipher_A0_I, decryptor, encoder,n);
    // Ciphertext cipher_Y7_tmp_change = pack_change_cipher_to_matrix_A_new_encoding(
    //     context,cipher_Y7, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // // decrypt_and_print_new_encoding(cipher_A_change, decryptor, encoder,n);
    // Ciphertext cipher_A7_change_I = pack_change_cipher_to_matrix_B_new_encoding(
    //     context,cipher_A7_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // // decrypt_and_print_new_encoding(cipher_B_change, decryptor, encoder,n);

    // Ciphertext cipher_Y8 = pack_encrypted_matrix_mul_new_encoding(
    //     context,cipher_Y7_tmp_change,cipher_A7_change_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);

    // Ciphertext cipher_A7_change_1 = pack_change_cipher_to_matrix_A_new_encoding(
    //     context,cipher_A7, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // Ciphertext cipher_A7_change_2 = pack_change_cipher_to_matrix_B_new_encoding(
    //     context,cipher_A7, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    // Ciphertext cipher_A8 = pack_encrypted_matrix_mul_new_encoding(
    //     context,cipher_A7_change_1,cipher_A7_change_2, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    
    // decrypt_and_print_new_encoding(cipher_Y8, decryptor, encoder,n);
    // auto end_time1 = chrono::high_resolution_clock::now();
    // auto time_diff1 = chrono::duration_cast<chrono::microseconds>(end_time1 - start_time1);


    // auto time_all = time_diff0+time_diff1;
    // cout <<"the all time is :"<<time_all.count() / 1e6 << " s" << endl;
    // std::cout << std::endl;

    // seal::Plaintext plain_I_8;
    // seal::Ciphertext cipher_A8_I;
    // encoder.encode(result_I,cipher_A8.scale(),plain_I_8);
    // evaluator.mod_switch_to_inplace(plain_I_8,cipher_A8.parms_id());
    // evaluator.add_plain(cipher_A8,plain_I_8,cipher_A8_I);

    // Ciphertext cipher_Y8_change = pack_change_cipher_to_matrix_A_new_encoding(
    //     context,cipher_Y8, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // // decrypt_and_print_new_encoding(cipher_A_change, decryptor, encoder,n);


    // Ciphertext cipher_A8_change_I = pack_change_cipher_to_matrix_B_new_encoding(
    //     context,cipher_A8_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // // decrypt_and_print_new_encoding(cipher_B_change, decryptor, encoder,n);
    // Ciphertext cipher_Y9 = pack_encrypted_matrix_mul_new_encoding(
    //     context,cipher_Y8_change,cipher_A8_change_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);

    // Ciphertext cipher_A8_change_1 = pack_change_cipher_to_matrix_A_new_encoding(
    //     context,cipher_A8, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // Ciphertext cipher_A8_change_2 = pack_change_cipher_to_matrix_B_new_encoding(
    //     context,cipher_A8, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    // Ciphertext cipher_A9 = pack_encrypted_matrix_mul_new_encoding(
    //     context,cipher_A8_change_1,cipher_A8_change_2, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    // // decrypt_and_print_new_encoding(cipher_Y9, decryptor, encoder,n);

    // auto end_time1 = chrono::high_resolution_clock::now();
    // auto time_diff1= chrono::duration_cast<chrono::microseconds>(end_time1- start_time1);



    // // // //到这里是r=9


    // // auto time_all = time_diff0+time_diff1;
    // // cout <<"the all time is :"<<time_all.count() / 1e6 << " s" << endl;

    // seal:: Plaintext plain_A9,plain_Y9,plain_A9_1,plain_Y9_1;
    // seal:: Ciphertext cipher_A9_tmp,cipher_Y9_tmp;
    // decryptor.decrypt(cipher_A9,plain_A9);
    // decryptor.decrypt(cipher_Y9,plain_Y9);

    // std::vector<double> result_A9;
    // std::vector<double> result_Y9;
    // encoder.decode(plain_A9, result_A9);
    // encoder.decode(plain_Y9, result_Y9);

    // auto start_time2 = chrono::high_resolution_clock::now();

    // encoder.encode(result_A9,scale,plain_A9_1);
    // encoder.encode(result_Y9,scale,plain_Y9_1);

    // encryptor.encrypt(plain_A9_1,cipher_A9_tmp);
    // encryptor.encrypt(plain_Y9_1,cipher_Y9_tmp);

    // seal::Plaintext plain_I_9;
    // seal::Ciphertext cipher_A9_I;
    // encoder.encode(result_I,cipher_A9_tmp.scale(),plain_I_9);
    // evaluator.mod_switch_to_inplace(plain_I_9,cipher_A9_tmp.parms_id());
    // evaluator.add_plain(cipher_A9_tmp,plain_I_9,cipher_A9_I);
    // // decrypt_and_print_new_encoding(cipher_A0_I, decryptor, encoder,n);
    // Ciphertext cipher_Y9_tmp_change = pack_change_cipher_to_matrix_A_new_encoding(
    //     context,cipher_Y9_tmp, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // // decrypt_and_print_new_encoding(cipher_A_change, decryptor, encoder,n);
    // Ciphertext cipher_A9_change_I = pack_change_cipher_to_matrix_B_new_encoding(
    //     context,cipher_A9_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // // decrypt_and_print_new_encoding(cipher_B_change, decryptor, encoder,n);

    // Ciphertext cipher_Y10 = pack_encrypted_matrix_mul_new_encoding(
    //     context,cipher_Y9_tmp_change,cipher_A9_change_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);

    // Ciphertext cipher_A9_change_1 = pack_change_cipher_to_matrix_A_new_encoding(
    //     context,cipher_A9_tmp, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // Ciphertext cipher_A9_change_2 = pack_change_cipher_to_matrix_B_new_encoding(
    //     context,cipher_A9_tmp, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    // Ciphertext cipher_A10 = pack_encrypted_matrix_mul_new_encoding(
    //     context,cipher_A9_change_1,cipher_A9_change_2, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    // // decrypt_and_print_new_encoding(cipher_Y10, decryptor, encoder,n);

    // // auto end_time2 = chrono::high_resolution_clock::now();
    // // auto time_diff2= chrono::duration_cast<chrono::microseconds>(end_time2- start_time2);



    // // // //到这里是r=10


    // // auto time_all = time_diff0+time_diff1+time_diff2;
    // // cout <<"the all time is :"<<time_all.count() / 1e6 << " s" << endl;
    // seal::Plaintext plain_I_10;
    // seal::Ciphertext cipher_A10_I;
    // encoder.encode(result_I,cipher_A10.scale(),plain_I_10);
    // evaluator.mod_switch_to_inplace(plain_I_10,cipher_A10.parms_id());
    // evaluator.add_plain(cipher_A10,plain_I_10,cipher_A10_I);

    // Ciphertext cipher_Y10_change = pack_change_cipher_to_matrix_A_new_encoding(
    //     context,cipher_Y10, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // // decrypt_and_print_new_encoding(cipher_A_change, decryptor, encoder,n);


    // Ciphertext cipher_A10_change_I = pack_change_cipher_to_matrix_B_new_encoding(
    //     context,cipher_A10_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // // decrypt_and_print_new_encoding(cipher_B_change, decryptor, encoder,n);
    // Ciphertext cipher_Y11 = pack_encrypted_matrix_mul_new_encoding(
    //     context,cipher_Y10_change,cipher_A10_change_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);

    // Ciphertext cipher_A10_change_1 = pack_change_cipher_to_matrix_A_new_encoding(
    //     context,cipher_A10, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // Ciphertext cipher_A10_change_2 = pack_change_cipher_to_matrix_B_new_encoding(
    //     context,cipher_A10, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    // Ciphertext cipher_A11 = pack_encrypted_matrix_mul_new_encoding(
    //     context,cipher_A10_change_1,cipher_A10_change_2, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    // decrypt_and_print_new_encoding(cipher_Y11, decryptor, encoder,n);
    // auto end_time2  = chrono::high_resolution_clock::now();
    // auto time_diff2= chrono::duration_cast<chrono::microseconds>(end_time2 - start_time2);
    // // //到这里是r=11
    // auto time_all = time_diff0+time_diff1+time_diff2;
    // cout <<"the all time is :"<<time_all.count() / 1e6 << " s" << endl;

    // seal::Plaintext plain_I_11;
    // seal::Ciphertext cipher_A11_I;
    // encoder.encode(result_I,cipher_A11.scale(),plain_I_11);
    // evaluator.mod_switch_to_inplace(plain_I_11,cipher_A11.parms_id());
    // evaluator.add_plain(cipher_A11,plain_I_11,cipher_A11_I);
    // // decrypt_and_print_new_encoding(cipher_A0_I, decryptor, encoder,n);
    // Ciphertext cipher_Y11_tmp_change = pack_change_cipher_to_matrix_A_new_encoding(
    //     context,cipher_Y11, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // // decrypt_and_print_new_encoding(cipher_A_change, decryptor, encoder,n);
    // Ciphertext cipher_A11_change_I = pack_change_cipher_to_matrix_B_new_encoding(
    //     context,cipher_A11_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // // decrypt_and_print_new_encoding(cipher_B_change, decryptor, encoder,n);

    // Ciphertext cipher_Y12 = pack_encrypted_matrix_mul_new_encoding(
    //     context,cipher_Y11_tmp_change,cipher_A11_change_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);

    // Ciphertext cipher_A11_change_1 = pack_change_cipher_to_matrix_A_new_encoding(
    //     context,cipher_A11, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // Ciphertext cipher_A11_change_2 = pack_change_cipher_to_matrix_B_new_encoding(
    //     context,cipher_A11, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    // Ciphertext cipher_A12 = pack_encrypted_matrix_mul_new_encoding(
    //     context,cipher_A11_change_1,cipher_A11_change_2, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    // seal::Plaintext plain_I_12;
    // seal::Ciphertext cipher_A12_I;
    // encoder.encode(result_I,cipher_A12.scale(),plain_I_12);
    // evaluator.mod_switch_to_inplace(plain_I_12,cipher_A12.parms_id());
    // evaluator.add_plain(cipher_A12,plain_I_12,cipher_A12_I);

    // Ciphertext cipher_Y12_change = pack_change_cipher_to_matrix_A_new_encoding(
    //     context,cipher_Y12, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // // decrypt_and_print_new_encoding(cipher_A_change, decryptor, encoder,n);


    // Ciphertext cipher_A12_change_I = pack_change_cipher_to_matrix_B_new_encoding(
    //     context,cipher_A12_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // // decrypt_and_print_new_encoding(cipher_B_change, decryptor, encoder,n);
    // Ciphertext cipher_Y13 = pack_encrypted_matrix_mul_new_encoding(
    //     context,cipher_Y12_change,cipher_A12_change_I, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);

    // Ciphertext cipher_A12_change_1 = pack_change_cipher_to_matrix_A_new_encoding(
    //     context,cipher_A12, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    // Ciphertext cipher_A12_change_2 = pack_change_cipher_to_matrix_B_new_encoding(
    //     context,cipher_A12, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    // Ciphertext cipher_A13 = pack_encrypted_matrix_mul_new_encoding(
    //     context,cipher_A12_change_1,cipher_A12_change_2, scale, encoder, evaluator, decryptor, gal_keys, relin_keys,n);
    
    // decrypt_and_print_new_encoding(cipher_Y13, decryptor, encoder,n);

    // auto end_time2  = chrono::high_resolution_clock::now();
    // auto time_diff2= chrono::duration_cast<chrono::microseconds>(end_time2 - start_time2);
   

    // //到这里是r=13

    // auto time_all = time_diff0+time_diff1+time_diff2;
    // cout <<"the all time is :"<<time_all.count() / 1e6 << " s" << endl;
    // std::cout << std::endl;

}

void run_dia_block_inversion_512(size_t n, size_t m, size_t p,size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
//     parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30,30,30, 60}));
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 40, 40, 40, 40, 40, 40, 40, 40, 60}));
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 45, 45, 45, 45, 45, 45, 45, 45, 60}));
    double scale = pow(2.0, 45);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    // print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);
    cout<<"---------"<<endl;
    test_dia_block_inversion_512(context,n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}
void test_dia_block_inversion_512(SEALContext context,size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    // Matrix<double>A(n,m);
    // random_block_matrix_generator_512(n, m, A);
    // Matrix<double> C;//A^2
    // C = A.multiply(A);

    // Matrix<double>D;//A^2+A
    // D = C.add(A);
    size_t d=sqrt(encoder.slot_count());

    // Matrix<double>E;//A^3+A^4
    // E = C.multiply(D);

    // Matrix<double>F;//A^2+A+A^2+A
    // F = E.add(D);

    // std::ofstream outFile("matrix_plain_result_512.txt");
    // for (size_t i = 0; i < n; i++) {
    //     for (size_t j = 0; j < m; j++) {
    //         outFile << F.get(j, i) << (j < n - 1 ? " " : "");
    //     }
    //     outFile << std::endl;
    // }

    // outFile.close();
    Matrix<double>A(n,m);
    std::string filename = "../build/dia_input/dia_input_c_(512,512).txt";
    A.readFromFile(filename);
    cout<<"   wssssssssss"<<endl;
    // cout<<A.get(0,1);

    A.print(); // 打印矩阵以验证

    Matrix<double>A_0(d,d);
    Matrix<double>A_1(d,d);
    Matrix<double>A_2(d,d);
    Matrix<double>A_3(d,d);
    Matrix<double>A_4(d,d);
    Matrix<double>A_5(d,d);
    Matrix<double>A_6(d,d);
    Matrix<double>A_7(d,d);
    Matrix<double>A_8(d,d);
    Matrix<double>A_9(d,d);
    Matrix<double>A_10(d,d);
    Matrix<double>A_11(d,d);
    Matrix<double>A_12(d,d);
    Matrix<double>A_13(d,d);
    Matrix<double>A_14(d,d);
    Matrix<double>A_15(d,d);

    A.split(d,A_0, A_1,A_2,A_3,A_4,A_5,A_6,A_7,A_8,A_9,A_10,A_11,A_12,A_13,A_14,A_15);


    auto start_time = chrono::high_resolution_clock::now();

    cipher_matrix_jiang cipher_A0(A_0,scale,encoder,encryptor);//A_0
    cipher_matrix_jiang cipher_A1(A_1,scale,encoder,encryptor);//A_1
    cipher_matrix_jiang cipher_A2(A_2,scale,encoder,encryptor);//A_2
    cipher_matrix_jiang cipher_A3(A_3,scale,encoder,encryptor);//A_3
    cipher_matrix_jiang cipher_A4(A_4,scale,encoder,encryptor);//A_0
    cipher_matrix_jiang cipher_A5(A_5,scale,encoder,encryptor);//A_1
    cipher_matrix_jiang cipher_A6(A_6,scale,encoder,encryptor);//A_2
    cipher_matrix_jiang cipher_A7(A_7,scale,encoder,encryptor);//A_3
    cipher_matrix_jiang cipher_A8(A_8,scale,encoder,encryptor);//A_0
    cipher_matrix_jiang cipher_A9(A_9,scale,encoder,encryptor);//A_1
    cipher_matrix_jiang cipher_A10(A_10,scale,encoder,encryptor);//A_2
    cipher_matrix_jiang cipher_A11(A_11,scale,encoder,encryptor);//A_3
    cipher_matrix_jiang cipher_A12(A_12,scale,encoder,encryptor);//A_0
    cipher_matrix_jiang cipher_A13(A_13,scale,encoder,encryptor);//A_1
    cipher_matrix_jiang cipher_A14(A_14,scale,encoder,encryptor);//A_2
    cipher_matrix_jiang cipher_A15(A_15,scale,encoder,encryptor);//A_3

    cipher_matrix_jiang cipher_B0;
    cipher_matrix_jiang cipher_B1;
    cipher_matrix_jiang cipher_B2;
    cipher_matrix_jiang cipher_B3;
    cipher_matrix_jiang cipher_B4;
    cipher_matrix_jiang cipher_B5;
    cipher_matrix_jiang cipher_B6;
    cipher_matrix_jiang cipher_B7;
    cipher_matrix_jiang cipher_B8;
    cipher_matrix_jiang cipher_B9;
    cipher_matrix_jiang cipher_B10;
    cipher_matrix_jiang cipher_B11;
    cipher_matrix_jiang cipher_B12;
    cipher_matrix_jiang cipher_B13;
    cipher_matrix_jiang cipher_B14;
    cipher_matrix_jiang cipher_B15;
    cipher_B0=cipher_A0;
    cipher_B1=cipher_A1;
    cipher_B2=cipher_A2;
    cipher_B3=cipher_A3;
    cipher_B4=cipher_A4;
    cipher_B5=cipher_A5;
    cipher_B6=cipher_A6;
    cipher_B7=cipher_A7;
    cipher_B8=cipher_A8;
    cipher_B9=cipher_A9;
    cipher_B10=cipher_A10;
    cipher_B11=cipher_A11;
    cipher_B12=cipher_A12;
    cipher_B13=cipher_A13;
    cipher_B14=cipher_A14;
    cipher_B15=cipher_A15;

    seal::Ciphertext cipher_A0_ciphertext,cipher_A1_ciphertext,cipher_A2_ciphertext,cipher_A3_ciphertext;
    seal::Ciphertext cipher_A4_ciphertext,cipher_A5_ciphertext,cipher_A6_ciphertext,cipher_A7_ciphertext;
    seal::Ciphertext cipher_A8_ciphertext,cipher_A9_ciphertext,cipher_A10_ciphertext,cipher_A11_ciphertext;
    seal::Ciphertext cipher_A12_ciphertext,cipher_A13_ciphertext,cipher_A14_ciphertext,cipher_A15_ciphertext;
 

    cipher_A0_ciphertext=cipher_A0.get_cipher_matrix();
    cipher_A1_ciphertext=cipher_A1.get_cipher_matrix();
    cipher_A2_ciphertext=cipher_A2.get_cipher_matrix();
    cipher_A3_ciphertext=cipher_A3.get_cipher_matrix();
    cipher_A4_ciphertext=cipher_A4.get_cipher_matrix();
    cipher_A5_ciphertext=cipher_A5.get_cipher_matrix();
    cipher_A6_ciphertext=cipher_A6.get_cipher_matrix();
    cipher_A7_ciphertext=cipher_A7.get_cipher_matrix();
    cipher_A8_ciphertext=cipher_A8.get_cipher_matrix();
    cipher_A9_ciphertext=cipher_A9.get_cipher_matrix();
    cipher_A10_ciphertext=cipher_A10.get_cipher_matrix();
    cipher_A11_ciphertext=cipher_A11.get_cipher_matrix();
    cipher_A12_ciphertext=cipher_A12.get_cipher_matrix();
    cipher_A13_ciphertext=cipher_A13.get_cipher_matrix();
    cipher_A14_ciphertext=cipher_A14.get_cipher_matrix();
    cipher_A15_ciphertext=cipher_A15.get_cipher_matrix();

    cipher_A0.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B0.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A1.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B1.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A2.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B2.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A3.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B3.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A4.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B4.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A5.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B5.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A6.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B6.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A7.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B7.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A8.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B8.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A9.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B9.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A10.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B10.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A11.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B11.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A12.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B12.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A13.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B13.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A14.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B14.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_A15.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_B15.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    
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

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A0,cipher_B0,cipher_A0B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A0,cipher_B1,cipher_A0B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A0,cipher_B2,cipher_A0B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A0,cipher_B3,cipher_A0B3);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A4,cipher_B0,cipher_A4B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A4,cipher_B1,cipher_A4B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A4,cipher_B2,cipher_A4B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A4,cipher_B3,cipher_A4B3);
        
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A8,cipher_B0,cipher_A8B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A8,cipher_B1,cipher_A8B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A8,cipher_B2,cipher_A8B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A8,cipher_B3,cipher_A8B3);
        
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A12,cipher_B0,cipher_A12B0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A12,cipher_B1,cipher_A12B1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A12,cipher_B2,cipher_A12B2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A12,cipher_B3,cipher_A12B3);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1,cipher_B4,cipher_A1B4);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1,cipher_B5,cipher_A1B5);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1,cipher_B6,cipher_A1B6);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A1,cipher_B7,cipher_A1B7);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A5,cipher_B4,cipher_A5B4);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A5,cipher_B5,cipher_A5B5);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A5,cipher_B6,cipher_A5B6);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A5,cipher_B7,cipher_A5B7);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A9,cipher_B4,cipher_A9B4);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A9,cipher_B5,cipher_A9B5);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A9,cipher_B6,cipher_A9B6);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A9,cipher_B7,cipher_A9B7);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A13,cipher_B4,cipher_A13B4);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A13,cipher_B5,cipher_A13B5);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A13,cipher_B6,cipher_A13B6);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A13,cipher_B7,cipher_A13B7);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A2,cipher_B8,cipher_A2B8);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A2,cipher_B9,cipher_A2B9);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A2,cipher_B10,cipher_A2B10);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A2,cipher_B11,cipher_A2B11);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A6,cipher_B8,cipher_A6B8);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A6,cipher_B9,cipher_A6B9);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A6,cipher_B10,cipher_A6B10);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A6,cipher_B11,cipher_A6B11);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A10,cipher_B8,cipher_A10B8);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A10,cipher_B9,cipher_A10B9);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A10,cipher_B10,cipher_A10B10);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A10,cipher_B11,cipher_A10B11);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A14,cipher_B8,cipher_A14B8);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A14,cipher_B9,cipher_A14B9);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A14,cipher_B10,cipher_A14B10);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A14,cipher_B11,cipher_A14B11);
       
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A3,cipher_B12,cipher_A3B12);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A3,cipher_B13,cipher_A3B13);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A3,cipher_B14,cipher_A3B14);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A3,cipher_B15,cipher_A3B15);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A7,cipher_B12,cipher_A7B12);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A7,cipher_B13,cipher_A7B13);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A7,cipher_B14,cipher_A7B14);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A7,cipher_B15,cipher_A7B15);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A11,cipher_B12,cipher_A11B12);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A11,cipher_B13,cipher_A11B13);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A11,cipher_B14,cipher_A11B14);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A11,cipher_B15,cipher_A11B15);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A15,cipher_B12,cipher_A15B12);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A15,cipher_B13,cipher_A15B13);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A15,cipher_B14,cipher_A15B14);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_A15,cipher_B15,cipher_A15B15);
   
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

    cipher_matrix_jiang cipher_C0(cipher_temp_0,d,d,d);
    cipher_matrix_jiang cipher_C1(cipher_temp_1,d,d,d);
    cipher_matrix_jiang cipher_C2(cipher_temp_2,d,d,d);
    cipher_matrix_jiang cipher_C3(cipher_temp_3,d,d,d);
    cipher_matrix_jiang cipher_C4(cipher_temp_4,d,d,d);
    cipher_matrix_jiang cipher_C5(cipher_temp_5,d,d,d);
    cipher_matrix_jiang cipher_C6(cipher_temp_6,d,d,d);
    cipher_matrix_jiang cipher_C7(cipher_temp_7,d,d,d);
    cipher_matrix_jiang cipher_C8(cipher_temp_8,d,d,d);
    cipher_matrix_jiang cipher_C9(cipher_temp_9,d,d,d);
    cipher_matrix_jiang cipher_C10(cipher_temp_10,d,d,d);
    cipher_matrix_jiang cipher_C11(cipher_temp_11,d,d,d);
    cipher_matrix_jiang cipher_C12(cipher_temp_12,d,d,d);
    cipher_matrix_jiang cipher_C13(cipher_temp_13,d,d,d);
    cipher_matrix_jiang cipher_C14(cipher_temp_14,d,d,d);
    cipher_matrix_jiang cipher_C15(cipher_temp_15,d,d,d);
    
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

    cipher_matrix_jiang cipher_D0(cipher_temp_0_A0,d,d,d);
    cipher_matrix_jiang cipher_D1(cipher_temp_1_A1,d,d,d);
    cipher_matrix_jiang cipher_D2(cipher_temp_2_A2,d,d,d);
    cipher_matrix_jiang cipher_D3(cipher_temp_3_A3,d,d,d);
    cipher_matrix_jiang cipher_D4(cipher_temp_4_A4,d,d,d);
    cipher_matrix_jiang cipher_D5(cipher_temp_5_A5,d,d,d);
    cipher_matrix_jiang cipher_D6(cipher_temp_6_A6,d,d,d);
    cipher_matrix_jiang cipher_D7(cipher_temp_7_A7,d,d,d);
    cipher_matrix_jiang cipher_D8(cipher_temp_8_A8,d,d,d);
    cipher_matrix_jiang cipher_D9(cipher_temp_9_A9,d,d,d);
    cipher_matrix_jiang cipher_D10(cipher_temp_10_A10,d,d,d);
    cipher_matrix_jiang cipher_D11(cipher_temp_11_A11,d,d,d);
    cipher_matrix_jiang cipher_D12(cipher_temp_12_A12,d,d,d);
    cipher_matrix_jiang cipher_D13(cipher_temp_13_A13,d,d,d);
    cipher_matrix_jiang cipher_D14(cipher_temp_14_A14,d,d,d);
    cipher_matrix_jiang cipher_D15(cipher_temp_15_A15,d,d,d);

    seal::Ciphertext cipher_C0_ciphertext,cipher_C1_ciphertext,cipher_C2_ciphertext,cipher_C3_ciphertext;
    seal::Ciphertext cipher_C4_ciphertext,cipher_C5_ciphertext,cipher_C6_ciphertext,cipher_C7_ciphertext;
    seal::Ciphertext cipher_C8_ciphertext,cipher_C9_ciphertext,cipher_C10_ciphertext,cipher_C11_ciphertext;
    seal::Ciphertext cipher_C12_ciphertext,cipher_C13_ciphertext,cipher_C14_ciphertext,cipher_C15_ciphertext;
 
    cipher_C0_ciphertext=cipher_C0.get_cipher_matrix();
    cipher_C1_ciphertext=cipher_C1.get_cipher_matrix();
    cipher_C2_ciphertext=cipher_C2.get_cipher_matrix();
    cipher_C3_ciphertext=cipher_C3.get_cipher_matrix();
    cipher_C4_ciphertext=cipher_C4.get_cipher_matrix();
    cipher_C5_ciphertext=cipher_C5.get_cipher_matrix();
    cipher_C6_ciphertext=cipher_C6.get_cipher_matrix();
    cipher_C7_ciphertext=cipher_C7.get_cipher_matrix();
    cipher_C8_ciphertext=cipher_C8.get_cipher_matrix();
    cipher_C9_ciphertext=cipher_C9.get_cipher_matrix();
    cipher_C10_ciphertext=cipher_C10.get_cipher_matrix();
    cipher_C11_ciphertext=cipher_C11.get_cipher_matrix();
    cipher_C12_ciphertext=cipher_C12.get_cipher_matrix();
    cipher_C13_ciphertext=cipher_C13.get_cipher_matrix();
    cipher_C14_ciphertext=cipher_C14.get_cipher_matrix();
    cipher_C15_ciphertext=cipher_C15.get_cipher_matrix();

    cipher_C0.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_D0.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_C1.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_D1.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_C2.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_D2.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_C3.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_D3.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_C4.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_D4.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_C5.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_D5.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_C6.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_D6.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_C7.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_D7.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_C8.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_D8.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_C9.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_D9.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_C10.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_D10.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_C11.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_D11.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_C12.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_D12.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_C13.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_D13.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_C14.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_D14.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_C15.change_cipher_to_matrix_A(scale,encoder,evaluator,gal_keys,relin_keys);
    cipher_D15.change_cipher_to_matrix_B(scale,encoder,evaluator,gal_keys,relin_keys);
    
    cipher_matrix_jiang cipher_C0D0,cipher_C1D4,cipher_C2D8,cipher_C3D12;
    cipher_matrix_jiang cipher_C0D1,cipher_C1D5,cipher_C2D9,cipher_C3D13;
    cipher_matrix_jiang cipher_C0D2,cipher_C1D6,cipher_C2D10,cipher_C3D14;
    cipher_matrix_jiang cipher_C0D3,cipher_C1D7,cipher_C2D11,cipher_C3D15;

    cipher_matrix_jiang cipher_C4D0,cipher_C5D4,cipher_C6D8,cipher_C7D12;
    cipher_matrix_jiang cipher_C4D1,cipher_C5D5,cipher_C6D9,cipher_C7D13;
    cipher_matrix_jiang cipher_C4D2,cipher_C5D6,cipher_C6D10,cipher_C7D14;
    cipher_matrix_jiang cipher_C4D3,cipher_C5D7,cipher_C6D11,cipher_C7D15;

    cipher_matrix_jiang cipher_C8D0,cipher_C9D4,cipher_C10D8,cipher_C11D12;
    cipher_matrix_jiang cipher_C8D1,cipher_C9D5,cipher_C10D9,cipher_C11D13;
    cipher_matrix_jiang cipher_C8D2,cipher_C9D6,cipher_C10D10,cipher_C11D14;
    cipher_matrix_jiang cipher_C8D3,cipher_C9D7,cipher_C10D11,cipher_C11D15;

    cipher_matrix_jiang cipher_C12D0,cipher_C13D4,cipher_C14D8,cipher_C15D12;
    cipher_matrix_jiang cipher_C12D1,cipher_C13D5,cipher_C14D9,cipher_C15D13;
    cipher_matrix_jiang cipher_C12D2,cipher_C13D6,cipher_C14D10,cipher_C15D14;
    cipher_matrix_jiang cipher_C12D3,cipher_C13D7,cipher_C14D11,cipher_C15D15;

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C0,cipher_D0,cipher_C0D0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C0,cipher_D1,cipher_C0D1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C0,cipher_D2,cipher_C0D2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C0,cipher_D3,cipher_C0D3);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C4,cipher_D0,cipher_C4D0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C4,cipher_D1,cipher_C4D1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C4,cipher_D2,cipher_C4D2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C4,cipher_D3,cipher_C4D3);
        
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C8,cipher_D0,cipher_C8D0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C8,cipher_D1,cipher_C8D1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C8,cipher_D2,cipher_C8D2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C8,cipher_D3,cipher_C8D3);
        
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C12,cipher_D0,cipher_C12D0);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C12,cipher_D1,cipher_C12D1);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C12,cipher_D2,cipher_C12D2);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C12,cipher_D3,cipher_C12D3);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C1,cipher_D4,cipher_C1D4);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C1,cipher_D5,cipher_C1D5);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C1,cipher_D6,cipher_C1D6);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C1,cipher_D7,cipher_C1D7);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C5,cipher_D4,cipher_C5D4);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C5,cipher_D5,cipher_C5D5);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C5,cipher_D6,cipher_C5D6);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C5,cipher_D7,cipher_C5D7);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C9,cipher_D4,cipher_C9D4);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C9,cipher_D5,cipher_C9D5);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C9,cipher_D6,cipher_C9D6);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C9,cipher_D7,cipher_C9D7);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C13,cipher_D4,cipher_C13D4);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C13,cipher_D5,cipher_C13D5);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C13,cipher_D6,cipher_C13D6);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C13,cipher_D7,cipher_C13D7);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C2,cipher_D8,cipher_C2D8);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C2,cipher_D9,cipher_C2D9);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C2,cipher_D10,cipher_C2D10);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C2,cipher_D11,cipher_C2D11);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C6,cipher_D8,cipher_C6D8);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C6,cipher_D9,cipher_C6D9);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C6,cipher_D10,cipher_C6D10);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C6,cipher_D11,cipher_C6D11);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C10,cipher_D8,cipher_C10D8);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C10,cipher_D9,cipher_C10D9);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C10,cipher_D10,cipher_C10D10);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C10,cipher_D11,cipher_C10D11);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C14,cipher_D8,cipher_C14D8);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C14,cipher_D9,cipher_C14D9);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C14,cipher_D10,cipher_C14D10);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C14,cipher_D11,cipher_C14D11);
       
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C3,cipher_D12,cipher_C3D12);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C3,cipher_D13,cipher_C3D13);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C3,cipher_D14,cipher_C3D14);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C3,cipher_D15,cipher_C3D15);
    
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C7,cipher_D12,cipher_C7D12);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C7,cipher_D13,cipher_C7D13);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C7,cipher_D14,cipher_C7D14);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C7,cipher_D15,cipher_C7D15);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C11,cipher_D12,cipher_C11D12);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C11,cipher_D13,cipher_C11D13);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C11,cipher_D14,cipher_C11D14);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C11,cipher_D15,cipher_C11D15);

    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C15,cipher_D12,cipher_C15D12);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C15,cipher_D13,cipher_C15D13);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C15,cipher_D14,cipher_C15D14);
    encrypted_jiang_matrix_multiplication(encoder,scale,evaluator,relin_keys,gal_keys,cipher_C15,cipher_D15,cipher_C15D15);
   
    seal::Ciphertext cipher_C0D0_ciphertext,cipher_C1D4_ciphertext,cipher_C2D8_ciphertext,cipher_C3D12_ciphertext;
    seal::Ciphertext cipher_C0D1_ciphertext,cipher_C1D5_ciphertext,cipher_C2D9_ciphertext,cipher_C3D13_ciphertext;
    seal::Ciphertext cipher_C0D2_ciphertext,cipher_C1D6_ciphertext,cipher_C2D10_ciphertext,cipher_C3D14_ciphertext;
    seal::Ciphertext cipher_C0D3_ciphertext,cipher_C1D7_ciphertext,cipher_C2D11_ciphertext,cipher_C3D15_ciphertext;
    seal::Ciphertext cipher_C4D0_ciphertext,cipher_C5D4_ciphertext,cipher_C6D8_ciphertext,cipher_C7D12_ciphertext;
    seal::Ciphertext cipher_C4D1_ciphertext,cipher_C5D5_ciphertext,cipher_C6D9_ciphertext,cipher_C7D13_ciphertext;
    seal::Ciphertext cipher_C4D2_ciphertext,cipher_C5D6_ciphertext,cipher_C6D10_ciphertext,cipher_C7D14_ciphertext;
    seal::Ciphertext cipher_C4D3_ciphertext,cipher_C5D7_ciphertext,cipher_C6D11_ciphertext,cipher_C7D15_ciphertext;

    seal::Ciphertext cipher_C8D0_ciphertext,cipher_C9D4_ciphertext,cipher_C10D8_ciphertext,cipher_C11D12_ciphertext;
    seal::Ciphertext cipher_C8D1_ciphertext,cipher_C9D5_ciphertext,cipher_C10D9_ciphertext,cipher_C11D13_ciphertext;
    seal::Ciphertext cipher_C8D2_ciphertext,cipher_C9D6_ciphertext,cipher_C10D10_ciphertext,cipher_C11D14_ciphertext;
    seal::Ciphertext cipher_C8D3_ciphertext,cipher_C9D7_ciphertext,cipher_C10D11_ciphertext,cipher_C11D15_ciphertext;

    seal::Ciphertext cipher_C12D0_ciphertext,cipher_C13D4_ciphertext,cipher_C14D8_ciphertext,cipher_C15D12_ciphertext;
    seal::Ciphertext cipher_C12D1_ciphertext,cipher_C13D5_ciphertext,cipher_C14D9_ciphertext,cipher_C15D13_ciphertext;
    seal::Ciphertext cipher_C12D2_ciphertext,cipher_C13D6_ciphertext,cipher_C14D10_ciphertext,cipher_C15D14_ciphertext;
    seal::Ciphertext cipher_C12D3_ciphertext,cipher_C13D7_ciphertext,cipher_C14D11_ciphertext,cipher_C15D15_ciphertext;
    
    seal::Ciphertext cipher_temp_1_0,cipher_temp_1_1,cipher_temp_1_2,cipher_temp_1_3;
    seal::Ciphertext cipher_temp_1_4,cipher_temp_1_5,cipher_temp_1_6,cipher_temp_1_7;
    seal::Ciphertext cipher_temp_1_8,cipher_temp_1_9,cipher_temp_1_10,cipher_temp_1_11;
    seal::Ciphertext cipher_temp_1_12,cipher_temp_1_13,cipher_temp_1_14,cipher_temp_1_15;

    cipher_C0D0_ciphertext=cipher_C0D0.get_cipher_matrix();
    cipher_C0D1_ciphertext=cipher_C0D1.get_cipher_matrix();
    cipher_C0D2_ciphertext=cipher_C0D2.get_cipher_matrix();
    cipher_C0D3_ciphertext=cipher_C0D3.get_cipher_matrix();
    cipher_C4D0_ciphertext=cipher_C4D0.get_cipher_matrix();
    cipher_C4D1_ciphertext=cipher_C4D1.get_cipher_matrix();
    cipher_C4D2_ciphertext=cipher_C4D2.get_cipher_matrix();
    cipher_C4D3_ciphertext=cipher_C4D3.get_cipher_matrix();
    cipher_C8D0_ciphertext=cipher_C8D0.get_cipher_matrix();
    cipher_C8D1_ciphertext=cipher_C8D1.get_cipher_matrix();
    cipher_C8D2_ciphertext=cipher_C8D2.get_cipher_matrix();
    cipher_C8D3_ciphertext=cipher_C8D3.get_cipher_matrix();
    cipher_C12D0_ciphertext=cipher_C12D0.get_cipher_matrix();
    cipher_C12D1_ciphertext=cipher_C12D1.get_cipher_matrix();
    cipher_C12D2_ciphertext=cipher_C12D2.get_cipher_matrix();
    cipher_C12D3_ciphertext=cipher_C12D3.get_cipher_matrix();

    cipher_C1D4_ciphertext=cipher_C1D4.get_cipher_matrix();
    cipher_C1D5_ciphertext=cipher_C1D5.get_cipher_matrix();
    cipher_C1D6_ciphertext=cipher_C1D6.get_cipher_matrix();
    cipher_C1D7_ciphertext=cipher_C1D7.get_cipher_matrix();
    cipher_C5D4_ciphertext=cipher_C5D4.get_cipher_matrix();
    cipher_C5D5_ciphertext=cipher_C5D5.get_cipher_matrix();
    cipher_C5D6_ciphertext=cipher_C5D6.get_cipher_matrix();
    cipher_C5D7_ciphertext=cipher_C5D7.get_cipher_matrix();
    cipher_C9D4_ciphertext=cipher_C9D4.get_cipher_matrix();
    cipher_C9D5_ciphertext=cipher_C9D5.get_cipher_matrix();
    cipher_C9D6_ciphertext=cipher_C9D6.get_cipher_matrix();
    cipher_C9D7_ciphertext=cipher_C9D7.get_cipher_matrix();
    cipher_C13D4_ciphertext=cipher_C13D4.get_cipher_matrix();
    cipher_C13D5_ciphertext=cipher_C13D5.get_cipher_matrix();
    cipher_C13D6_ciphertext=cipher_C13D6.get_cipher_matrix();
    cipher_C13D7_ciphertext=cipher_C13D7.get_cipher_matrix();

    cipher_C2D8_ciphertext=cipher_C2D8.get_cipher_matrix();
    cipher_C2D9_ciphertext=cipher_C2D9.get_cipher_matrix();
    cipher_C2D10_ciphertext=cipher_C2D10.get_cipher_matrix();
    cipher_C2D11_ciphertext=cipher_C2D11.get_cipher_matrix();
    cipher_C6D8_ciphertext=cipher_C6D8.get_cipher_matrix();
    cipher_C6D9_ciphertext=cipher_C6D9.get_cipher_matrix();
    cipher_C6D10_ciphertext=cipher_C6D10.get_cipher_matrix();
    cipher_C6D11_ciphertext=cipher_C6D11.get_cipher_matrix();
    cipher_C10D8_ciphertext=cipher_C10D8.get_cipher_matrix();
    cipher_C10D9_ciphertext=cipher_C10D9.get_cipher_matrix();
    cipher_C10D10_ciphertext=cipher_C10D10.get_cipher_matrix();
    cipher_C10D11_ciphertext=cipher_C10D11.get_cipher_matrix();
    cipher_C14D8_ciphertext=cipher_C14D8.get_cipher_matrix();
    cipher_C14D9_ciphertext=cipher_C14D9.get_cipher_matrix();
    cipher_C14D10_ciphertext=cipher_C14D10.get_cipher_matrix();
    cipher_C14D11_ciphertext=cipher_C14D11.get_cipher_matrix();

    cipher_C3D12_ciphertext=cipher_C3D12.get_cipher_matrix();
    cipher_C3D13_ciphertext=cipher_C3D13.get_cipher_matrix();
    cipher_C3D14_ciphertext=cipher_C3D14.get_cipher_matrix();
    cipher_C3D15_ciphertext=cipher_C3D15.get_cipher_matrix();
    cipher_C7D12_ciphertext=cipher_C7D12.get_cipher_matrix();
    cipher_C7D13_ciphertext=cipher_C7D13.get_cipher_matrix();
    cipher_C7D14_ciphertext=cipher_C7D14.get_cipher_matrix();
    cipher_C7D15_ciphertext=cipher_C7D15.get_cipher_matrix();
    cipher_C11D12_ciphertext=cipher_C11D12.get_cipher_matrix();
    cipher_C11D13_ciphertext=cipher_C11D13.get_cipher_matrix();
    cipher_C11D14_ciphertext=cipher_C11D14.get_cipher_matrix();
    cipher_C11D15_ciphertext=cipher_C11D15.get_cipher_matrix();
    cipher_C15D12_ciphertext=cipher_C15D12.get_cipher_matrix();
    cipher_C15D13_ciphertext=cipher_C15D13.get_cipher_matrix();
    cipher_C15D14_ciphertext=cipher_C15D14.get_cipher_matrix();
    cipher_C15D15_ciphertext=cipher_C15D15.get_cipher_matrix();

    evaluator.add(cipher_C0D0_ciphertext,cipher_C1D4_ciphertext,cipher_temp_1_0);
    evaluator.add_inplace(cipher_temp_1_0,cipher_C2D8_ciphertext);
    evaluator.add_inplace(cipher_temp_1_0,cipher_C3D12_ciphertext);
    evaluator.add(cipher_C0D1_ciphertext,cipher_C1D5_ciphertext,cipher_temp_1_1);
    evaluator.add_inplace(cipher_temp_1_1,cipher_C2D9_ciphertext);
    evaluator.add_inplace(cipher_temp_1_1,cipher_C3D13_ciphertext);
    evaluator.add(cipher_C0D2_ciphertext,cipher_C1D6_ciphertext,cipher_temp_1_2);
    evaluator.add_inplace(cipher_temp_1_2,cipher_C2D10_ciphertext);
    evaluator.add_inplace(cipher_temp_1_2,cipher_C3D14_ciphertext);
    evaluator.add(cipher_C0D3_ciphertext,cipher_C1D7_ciphertext,cipher_temp_1_3);
    evaluator.add_inplace(cipher_temp_1_3,cipher_C2D11_ciphertext);
    evaluator.add_inplace(cipher_temp_1_3,cipher_C3D15_ciphertext);

    evaluator.add(cipher_C4D0_ciphertext,cipher_C5D4_ciphertext,cipher_temp_1_4);
    evaluator.add_inplace(cipher_temp_1_4,cipher_C6D8_ciphertext);
    evaluator.add_inplace(cipher_temp_1_4,cipher_C7D12_ciphertext);
    evaluator.add(cipher_C4D1_ciphertext,cipher_C5D5_ciphertext,cipher_temp_1_5);
    evaluator.add_inplace(cipher_temp_1_5,cipher_C6D9_ciphertext);
    evaluator.add_inplace(cipher_temp_1_5,cipher_C7D13_ciphertext);
    evaluator.add(cipher_C4D2_ciphertext,cipher_C5D6_ciphertext,cipher_temp_1_6);
    evaluator.add_inplace(cipher_temp_1_6,cipher_C6D10_ciphertext);
    evaluator.add_inplace(cipher_temp_1_6,cipher_C7D14_ciphertext);
    evaluator.add(cipher_C4D3_ciphertext,cipher_C5D7_ciphertext,cipher_temp_1_7);
    evaluator.add_inplace(cipher_temp_1_7,cipher_C6D11_ciphertext);
    evaluator.add_inplace(cipher_temp_1_7,cipher_C7D15_ciphertext);

    evaluator.add(cipher_C8D0_ciphertext,cipher_C9D4_ciphertext,cipher_temp_1_8);
    evaluator.add_inplace(cipher_temp_1_8,cipher_C10D8_ciphertext);
    evaluator.add_inplace(cipher_temp_1_8,cipher_C11D12_ciphertext);
    evaluator.add(cipher_C8D1_ciphertext,cipher_C9D5_ciphertext,cipher_temp_1_9);
    evaluator.add_inplace(cipher_temp_1_9,cipher_C10D9_ciphertext);
    evaluator.add_inplace(cipher_temp_1_9,cipher_C11D13_ciphertext);
    evaluator.add(cipher_C8D2_ciphertext,cipher_C9D6_ciphertext,cipher_temp_1_10);
    evaluator.add_inplace(cipher_temp_1_10,cipher_C10D10_ciphertext);
    evaluator.add_inplace(cipher_temp_1_10,cipher_C11D14_ciphertext);
    evaluator.add(cipher_C8D3_ciphertext,cipher_C9D7_ciphertext,cipher_temp_1_11);
    evaluator.add_inplace(cipher_temp_1_11,cipher_C10D11_ciphertext);
    evaluator.add_inplace(cipher_temp_1_11,cipher_C11D15_ciphertext);

    evaluator.add(cipher_C12D0_ciphertext,cipher_C13D4_ciphertext,cipher_temp_1_12);
    evaluator.add_inplace(cipher_temp_1_12,cipher_C14D8_ciphertext);
    evaluator.add_inplace(cipher_temp_1_12,cipher_C15D12_ciphertext);
    evaluator.add(cipher_C12D1_ciphertext,cipher_C13D5_ciphertext,cipher_temp_1_13);
    evaluator.add_inplace(cipher_temp_1_13,cipher_C14D9_ciphertext);
    evaluator.add_inplace(cipher_temp_1_13,cipher_C15D13_ciphertext);
    evaluator.add(cipher_C12D2_ciphertext,cipher_C13D6_ciphertext,cipher_temp_1_14);
    evaluator.add_inplace(cipher_temp_1_14,cipher_C14D10_ciphertext);
    evaluator.add_inplace(cipher_temp_1_14,cipher_C15D14_ciphertext);
    evaluator.add(cipher_C12D3_ciphertext,cipher_C13D7_ciphertext,cipher_temp_1_15);
    evaluator.add_inplace(cipher_temp_1_15,cipher_C14D11_ciphertext);
    evaluator.add_inplace(cipher_temp_1_15,cipher_C15D15_ciphertext);

    seal::Ciphertext result0,result1,result2,result3;
    seal::Ciphertext result4,result5,result6,result7;
    seal::Ciphertext result8,result9,result10,result11;
    seal::Ciphertext result12,result13,result14,result15;

    //compute (A^2+A)  with (A^3+A^4)
    match_levels(context,cipher_temp_0_A0,cipher_temp_1_0,evaluator);
    evaluator.add(cipher_temp_0_A0,cipher_temp_1_0,result0);
    match_levels(context,cipher_temp_1_A1,cipher_temp_1_1,evaluator);
    evaluator.add(cipher_temp_1_A1,cipher_temp_1_1,result1);
    match_levels(context,cipher_temp_2_A2,cipher_temp_1_2,evaluator);
    evaluator.add(cipher_temp_2_A2,cipher_temp_1_2,result2);
    match_levels(context,cipher_temp_3_A3,cipher_temp_1_3,evaluator);
    evaluator.add(cipher_temp_3_A3,cipher_temp_1_3,result3);
        
    match_levels(context,cipher_temp_4_A4,cipher_temp_1_4,evaluator);
    evaluator.add(cipher_temp_4_A4,cipher_temp_1_4,result4);
    match_levels(context,cipher_temp_5_A5,cipher_temp_1_5,evaluator);
    evaluator.add(cipher_temp_5_A5,cipher_temp_1_5,result5);
    match_levels(context,cipher_temp_6_A6,cipher_temp_1_6,evaluator);
    evaluator.add(cipher_temp_6_A6,cipher_temp_1_6,result6);
    match_levels(context,cipher_temp_7_A7,cipher_temp_1_7,evaluator);
    evaluator.add(cipher_temp_7_A7,cipher_temp_1_7,result7);
    match_levels(context,cipher_temp_8_A8,cipher_temp_1_8,evaluator);
    evaluator.add(cipher_temp_8_A8,cipher_temp_1_8,result8);
    match_levels(context,cipher_temp_9_A9,cipher_temp_1_9,evaluator);
    evaluator.add(cipher_temp_9_A9,cipher_temp_1_9,result9);
    match_levels(context,cipher_temp_10_A10,cipher_temp_1_10,evaluator);
    evaluator.add(cipher_temp_10_A10,cipher_temp_1_10,result10);
    match_levels(context,cipher_temp_11_A11,cipher_temp_1_11,evaluator);
    evaluator.add(cipher_temp_11_A11,cipher_temp_1_11,result11);
    match_levels(context,cipher_temp_12_A12,cipher_temp_1_12,evaluator);
    evaluator.add(cipher_temp_12_A12,cipher_temp_1_12,result12);
    match_levels(context,cipher_temp_13_A13,cipher_temp_1_13,evaluator);
    evaluator.add(cipher_temp_13_A13,cipher_temp_1_13,result13);
    match_levels(context,cipher_temp_14_A14,cipher_temp_1_14,evaluator);
    evaluator.add(cipher_temp_14_A14,cipher_temp_1_14,result14);
    match_levels(context,cipher_temp_15_A15,cipher_temp_1_15,evaluator);
    evaluator.add(cipher_temp_15_A15,cipher_temp_1_15,result15);

    cipher_matrix_jiang result_0(result0,d,d,d);
    cipher_matrix_jiang result_1(result1,d,d,d);
    cipher_matrix_jiang result_2(result2,d,d,d);
    cipher_matrix_jiang result_3(result3,d,d,d);
    cipher_matrix_jiang result_4(result4,d,d,d);
    cipher_matrix_jiang result_5(result5,d,d,d);
    cipher_matrix_jiang result_6(result6,d,d,d);
    cipher_matrix_jiang result_7(result7,d,d,d);
    cipher_matrix_jiang result_8(result8,d,d,d);
    cipher_matrix_jiang result_9(result9,d,d,d);
    cipher_matrix_jiang result_10(result10,d,d,d);
    cipher_matrix_jiang result_11(result11,d,d,d);
    cipher_matrix_jiang result_12(result12,d,d,d);
    cipher_matrix_jiang result_13(result13,d,d,d);
    cipher_matrix_jiang result_14(result14,d,d,d);
    cipher_matrix_jiang result_15(result15,d,d,d);

    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;

    Matrix<double> computed0,computed1,computed2,computed3,
    computed4,computed5,computed6,computed7,
    computed8,computed9,computed10,computed11,
    computed12,computed13,computed14,computed15;
    result_0.dec_matrix_cipher(computed0,encoder,decryptor);
    result_1.dec_matrix_cipher(computed1,encoder,decryptor);
    result_2.dec_matrix_cipher(computed2,encoder,decryptor);
    result_3.dec_matrix_cipher(computed3,encoder,decryptor);
    result_4.dec_matrix_cipher(computed4,encoder,decryptor);
    result_5.dec_matrix_cipher(computed5,encoder,decryptor);
    result_6.dec_matrix_cipher(computed6,encoder,decryptor);
    result_7.dec_matrix_cipher(computed7,encoder,decryptor);
    result_8.dec_matrix_cipher(computed8,encoder,decryptor);
    result_9.dec_matrix_cipher(computed9,encoder,decryptor);
    result_10.dec_matrix_cipher(computed10,encoder,decryptor);
    result_11.dec_matrix_cipher(computed11,encoder,decryptor);
    result_12.dec_matrix_cipher(computed12,encoder,decryptor);
    result_13.dec_matrix_cipher(computed13,encoder,decryptor);
    result_14.dec_matrix_cipher(computed14,encoder,decryptor);
    result_15.dec_matrix_cipher(computed15,encoder,decryptor);
    // cout << "the computed C = A^2+A: " << endl;
    //  store_matrix(n,m,computed_C);
    // computed.print();
    std::ofstream outFile0("matrix_result_512_0.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile0 << computed0.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile0 << std::endl;
    }

    outFile0.close();
    std::ofstream outFile1("matrix_result_512_1.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile1 << computed1.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile1 << std::endl;
    }

    outFile1.close();
    std::ofstream outFile2("matrix_result_512_2.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile2 << computed2.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile2 << std::endl;
    }

    outFile2.close();
    std::ofstream outFile3("matrix_result_512_3.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile3 << computed3.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile3 << std::endl;
    }

    outFile3.close();
    std::ofstream outFile4("matrix_result_512_4.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile4 << computed4.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile4 << std::endl;
    }

    outFile4.close();
    std::ofstream outFile5("matrix_result_512_5.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile5 << computed5.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile5 << std::endl;
    }

    outFile5.close();
    std::ofstream outFile6("matrix_result_512_6.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile6 << computed6.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile6 << std::endl;
    }

    outFile6.close();
    std::ofstream outFile7("matrix_result_512_7.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile7 << computed7.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile7 << std::endl;
    }

    outFile7.close();
    std::ofstream outFile8("matrix_result_512_8.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile8 << computed8.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile8 << std::endl;
    }

    outFile8.close();
    std::ofstream outFile9("matrix_result_512_9.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile9 << computed9.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile9 << std::endl;
    }

    outFile9.close();
    std::ofstream outFile10("matrix_result_512_10.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile10 << computed10.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile10 << std::endl;
    }

    outFile10.close();
    std::ofstream outFile11("matrix_result_512_11.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile11 << computed11.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile11 << std::endl;
    }

    outFile11.close();
    std::ofstream outFile12("matrix_result_512_12.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile12 << computed12.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile12 << std::endl;
    }

    outFile12.close();
    std::ofstream outFile13("matrix_result_512_13.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile13 << computed13.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile13 << std::endl;
    }

    outFile13.close();
    std::ofstream outFile14("matrix_result_512_14.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile14 << computed14.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile14 << std::endl;
    }

    outFile14.close();
    std::ofstream outFile15("matrix_result_512_15.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile15 << computed15.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile15 << std::endl;
    }

    outFile15.close();


    cout<<"endl"<<"---------------------------------------------------------"<<endl;
    // cout << "the true C =  A^2+A: " << endl;
    // F.print();

}
void run_real_block_inversion_512(size_t n, size_t m, size_t p,size_t poly_modulus_degree)
{
    EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
//     parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 30,30,30, 60}));
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 40, 40, 40, 40, 40, 40, 40, 40, 60}));
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, {50, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 60}));
    double scale = pow(2.0, 45);
    // cout << "scale = " << scale << endl;

    SEALContext context(parms);
    // print_parameters(context);

    KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    // size_t slot_count = encoder.slot_count();
    // cout << "Number of slots: " << slot_count << endl;

    GaloisKeys gal_keys;
    keygen.create_galois_keys(gal_keys);
    cout<<"---------"<<endl;
    test_real_block_inversion_512(context,n, m, p, encoder, encryptor, scale, evaluator, relin_keys, gal_keys, decryptor);
}
void test_real_block_inversion_512(SEALContext context,size_t n, size_t m, size_t p, CKKSEncoder &encoder, Encryptor &encryptor, double scale, Evaluator &evaluator, RelinKeys &relin_keys, GaloisKeys &gal_keys, Decryptor &decryptor)
{
    // Matrix<double>A(n,m);
    // random_block_matrix_generator_512(n, m, A);
    // Matrix<double> C;//A^2
    // C = A.multiply(A);

    // Matrix<double>D;//A^2+A
    // D = C.add(A);
    size_t d=sqrt(encoder.slot_count());

    // Matrix<double>E;//A^3+A^4
    // E = C.multiply(D);

    // Matrix<double>F;//A^2+A+A^2+A
    // F = E.add(D);

    // std::ofstream outFile("matrix_plain_result_512.txt");
    // for (size_t i = 0; i < n; i++) {
    //     for (size_t j = 0; j < m; j++) {
    //         outFile << F.get(j, i) << (j < n - 1 ? " " : "");
    //     }
    //     outFile << std::endl;
    // }

    // outFile.close();
    Matrix<double>A(n,m);
    std::string filename = "../build/real_input/real_input_c_(512,512).txt";
    A.readFromFile(filename);
    cout<<"   wssssssssss"<<endl;
    // cout<<A.get(0,1);

    A.print(); // 打印矩阵以验证

    Matrix<double>A_0(d,d);
    Matrix<double>A_1(d,d);
    Matrix<double>A_2(d,d);
    Matrix<double>A_3(d,d);
    Matrix<double>A_4(d,d);
    Matrix<double>A_5(d,d);
    Matrix<double>A_6(d,d);
    Matrix<double>A_7(d,d);
    Matrix<double>A_8(d,d);
    Matrix<double>A_9(d,d);
    Matrix<double>A_10(d,d);
    Matrix<double>A_11(d,d);
    Matrix<double>A_12(d,d);
    Matrix<double>A_13(d,d);
    Matrix<double>A_14(d,d);
    Matrix<double>A_15(d,d);

    A.split(d,A_0, A_1,A_2,A_3,A_4,A_5,A_6,A_7,A_8,A_9,A_10,A_11,A_12,A_13,A_14,A_15);


    auto start_time = chrono::high_resolution_clock::now();

    cipher_matrix_jiang cipher_A0(A_0,scale,encoder,encryptor);//A_0
    cipher_matrix_jiang cipher_A1(A_1,scale,encoder,encryptor);//A_1
    cipher_matrix_jiang cipher_A2(A_2,scale,encoder,encryptor);//A_2
    cipher_matrix_jiang cipher_A3(A_3,scale,encoder,encryptor);//A_3
    cipher_matrix_jiang cipher_A4(A_4,scale,encoder,encryptor);//A_0
    cipher_matrix_jiang cipher_A5(A_5,scale,encoder,encryptor);//A_1
    cipher_matrix_jiang cipher_A6(A_6,scale,encoder,encryptor);//A_2
    cipher_matrix_jiang cipher_A7(A_7,scale,encoder,encryptor);//A_3
    cipher_matrix_jiang cipher_A8(A_8,scale,encoder,encryptor);//A_0
    cipher_matrix_jiang cipher_A9(A_9,scale,encoder,encryptor);//A_1
    cipher_matrix_jiang cipher_A10(A_10,scale,encoder,encryptor);//A_2
    cipher_matrix_jiang cipher_A11(A_11,scale,encoder,encryptor);//A_3
    cipher_matrix_jiang cipher_A12(A_12,scale,encoder,encryptor);//A_0
    cipher_matrix_jiang cipher_A13(A_13,scale,encoder,encryptor);//A_1
    cipher_matrix_jiang cipher_A14(A_14,scale,encoder,encryptor);//A_2
    cipher_matrix_jiang cipher_A15(A_15,scale,encoder,encryptor);//A_3

    cipher_matrix_jiang cipher_B0;
    cipher_matrix_jiang cipher_B1;
    cipher_matrix_jiang cipher_B2;
    cipher_matrix_jiang cipher_B3;
    cipher_matrix_jiang cipher_B4;
    cipher_matrix_jiang cipher_B5;
    cipher_matrix_jiang cipher_B6;
    cipher_matrix_jiang cipher_B7;
    cipher_matrix_jiang cipher_B8;
    cipher_matrix_jiang cipher_B9;
    cipher_matrix_jiang cipher_B10;
    cipher_matrix_jiang cipher_B11;
    cipher_matrix_jiang cipher_B12;
    cipher_matrix_jiang cipher_B13;
    cipher_matrix_jiang cipher_B14;
    cipher_matrix_jiang cipher_B15;
    cipher_B0=cipher_A0;
    cipher_B1=cipher_A1;
    cipher_B2=cipher_A2;
    cipher_B3=cipher_A3;
    cipher_B4=cipher_A4;
    cipher_B5=cipher_A5;
    cipher_B6=cipher_A6;
    cipher_B7=cipher_A7;
    cipher_B8=cipher_A8;
    cipher_B9=cipher_A9;
    cipher_B10=cipher_A10;
    cipher_B11=cipher_A11;
    cipher_B12=cipher_A12;
    cipher_B13=cipher_A13;
    cipher_B14=cipher_A14;
    cipher_B15=cipher_A15;
    // A^2
    vector<Ciphertext> left_degree_2=block_512_LeftSquare(context,encoder,scale,evaluator,relin_keys,gal_keys,cipher_A0,cipher_A1,cipher_A2,cipher_A3,
        cipher_A4,cipher_A5,cipher_A6,cipher_A7,cipher_A8,cipher_A9,cipher_A10,cipher_A11,
        cipher_A12,cipher_A13,cipher_A14,cipher_A15,cipher_B0,cipher_B1,cipher_B2,cipher_B3,
        cipher_B4,cipher_B5,cipher_B6,cipher_B7,cipher_B8,cipher_B9,cipher_B10,cipher_B11,
        cipher_B12,cipher_B13,cipher_B14,cipher_B15);
    //A+A*A
    vector<Ciphertext> right_degree_2=block_512_LeftMulRight_AddRight(context,encoder,scale,evaluator,relin_keys,gal_keys,
        cipher_A0,cipher_A1,cipher_A2,cipher_A3,
        cipher_A4,cipher_A5,cipher_A6,cipher_A7,cipher_A8,cipher_A9,cipher_A10,cipher_A11,
        cipher_A12,cipher_A13,cipher_A14,cipher_A15,cipher_B0,cipher_B1,cipher_B2,cipher_B3,
        cipher_B4,cipher_B5,cipher_B6,cipher_B7,cipher_B8,cipher_B9,cipher_B10,cipher_B11,
        cipher_B12,cipher_B13,cipher_B14,cipher_B15);

    cipher_matrix_jiang left_degree_2_0(left_degree_2[0],d,d,d);
    cipher_matrix_jiang left_degree_2_1(left_degree_2[1],d,d,d);
    cipher_matrix_jiang left_degree_2_2(left_degree_2[2],d,d,d);
    cipher_matrix_jiang left_degree_2_3(left_degree_2[3],d,d,d);
    cipher_matrix_jiang left_degree_2_4(left_degree_2[4],d,d,d);
    cipher_matrix_jiang left_degree_2_5(left_degree_2[5],d,d,d);
    cipher_matrix_jiang left_degree_2_6(left_degree_2[6],d,d,d);
    cipher_matrix_jiang left_degree_2_7(left_degree_2[7],d,d,d);
    cipher_matrix_jiang left_degree_2_8(left_degree_2[8],d,d,d);
    cipher_matrix_jiang left_degree_2_9(left_degree_2[9],d,d,d);
    cipher_matrix_jiang left_degree_2_10(left_degree_2[10],d,d,d);
    cipher_matrix_jiang left_degree_2_11(left_degree_2[11],d,d,d);
    cipher_matrix_jiang left_degree_2_12(left_degree_2[12],d,d,d);
    cipher_matrix_jiang left_degree_2_13(left_degree_2[13],d,d,d);
    cipher_matrix_jiang left_degree_2_14(left_degree_2[14],d,d,d);
    cipher_matrix_jiang left_degree_2_15(left_degree_2[15],d,d,d);

    cipher_matrix_jiang right_degree_2_0(right_degree_2[0],d,d,d);
    cipher_matrix_jiang right_degree_2_1(right_degree_2[1],d,d,d);
    cipher_matrix_jiang right_degree_2_2(right_degree_2[2],d,d,d);
    cipher_matrix_jiang right_degree_2_3(right_degree_2[3],d,d,d);
    cipher_matrix_jiang right_degree_2_4(right_degree_2[4],d,d,d);
    cipher_matrix_jiang right_degree_2_5(right_degree_2[5],d,d,d);
    cipher_matrix_jiang right_degree_2_6(right_degree_2[6],d,d,d);
    cipher_matrix_jiang right_degree_2_7(right_degree_2[7],d,d,d);
    cipher_matrix_jiang right_degree_2_8(right_degree_2[8],d,d,d);
    cipher_matrix_jiang right_degree_2_9(right_degree_2[9],d,d,d);
    cipher_matrix_jiang right_degree_2_10(right_degree_2[10],d,d,d);
    cipher_matrix_jiang right_degree_2_11(right_degree_2[11],d,d,d);
    cipher_matrix_jiang right_degree_2_12(right_degree_2[12],d,d,d);
    cipher_matrix_jiang right_degree_2_13(right_degree_2[13],d,d,d);
    cipher_matrix_jiang right_degree_2_14(right_degree_2[14],d,d,d);
    cipher_matrix_jiang right_degree_2_15(right_degree_2[15],d,d,d);

    cout<<" 2222222222222222222222222222222 "<<endl;

    // A^4
    vector<Ciphertext> left_degree_4=block_512_LeftSquare(context,encoder,scale,evaluator,relin_keys,gal_keys,
        left_degree_2_0,left_degree_2_1,left_degree_2_2,left_degree_2_3,left_degree_2_4,left_degree_2_5,
        left_degree_2_6,left_degree_2_7,left_degree_2_8,left_degree_2_9,left_degree_2_10,left_degree_2_11,
        left_degree_2_12,left_degree_2_13,left_degree_2_14,left_degree_2_15,left_degree_2_0,left_degree_2_1,left_degree_2_2,left_degree_2_3,left_degree_2_4,left_degree_2_5,
        left_degree_2_6,left_degree_2_7,left_degree_2_8,left_degree_2_9,left_degree_2_10,left_degree_2_11,
        left_degree_2_12,left_degree_2_13,left_degree_2_14,left_degree_2_15);
    //A+...+A^4
    vector<Ciphertext> right_degree_4=block_512_LeftMulRight_AddRight(context,encoder,scale,evaluator,relin_keys,gal_keys,
        left_degree_2_0,left_degree_2_1,left_degree_2_2,left_degree_2_3,left_degree_2_4,left_degree_2_5,
        left_degree_2_6,left_degree_2_7,left_degree_2_8,left_degree_2_9,left_degree_2_10,left_degree_2_11,
        left_degree_2_12,left_degree_2_13,left_degree_2_14,left_degree_2_15,right_degree_2_0,right_degree_2_1,right_degree_2_2,right_degree_2_3,right_degree_2_4,right_degree_2_5,
        right_degree_2_6,right_degree_2_7,right_degree_2_8,right_degree_2_9,right_degree_2_10,right_degree_2_11,
        right_degree_2_12,right_degree_2_13,right_degree_2_14,right_degree_2_15);
    
    cipher_matrix_jiang left_degree_4_0(left_degree_4[0],d,d,d);
    cipher_matrix_jiang left_degree_4_1(left_degree_4[1],d,d,d);
    cipher_matrix_jiang left_degree_4_2(left_degree_4[2],d,d,d);
    cipher_matrix_jiang left_degree_4_3(left_degree_4[3],d,d,d);
    cipher_matrix_jiang left_degree_4_4(left_degree_4[4],d,d,d);
    cipher_matrix_jiang left_degree_4_5(left_degree_4[5],d,d,d);
    cipher_matrix_jiang left_degree_4_6(left_degree_4[6],d,d,d);
    cipher_matrix_jiang left_degree_4_7(left_degree_4[7],d,d,d);
    cipher_matrix_jiang left_degree_4_8(left_degree_4[8],d,d,d);
    cipher_matrix_jiang left_degree_4_9(left_degree_4[9],d,d,d);
    cipher_matrix_jiang left_degree_4_10(left_degree_4[10],d,d,d);
    cipher_matrix_jiang left_degree_4_11(left_degree_4[11],d,d,d);
    cipher_matrix_jiang left_degree_4_12(left_degree_4[12],d,d,d);
    cipher_matrix_jiang left_degree_4_13(left_degree_4[13],d,d,d);
    cipher_matrix_jiang left_degree_4_14(left_degree_4[14],d,d,d);
    cipher_matrix_jiang left_degree_4_15(left_degree_4[15],d,d,d);

    cipher_matrix_jiang right_degree_4_0(right_degree_4[0],d,d,d);
    cipher_matrix_jiang right_degree_4_1(right_degree_4[1],d,d,d);
    cipher_matrix_jiang right_degree_4_2(right_degree_4[2],d,d,d);
    cipher_matrix_jiang right_degree_4_3(right_degree_4[3],d,d,d);
    cipher_matrix_jiang right_degree_4_4(right_degree_4[4],d,d,d);
    cipher_matrix_jiang right_degree_4_5(right_degree_4[5],d,d,d);
    cipher_matrix_jiang right_degree_4_6(right_degree_4[6],d,d,d);
    cipher_matrix_jiang right_degree_4_7(right_degree_4[7],d,d,d);
    cipher_matrix_jiang right_degree_4_8(right_degree_4[8],d,d,d);
    cipher_matrix_jiang right_degree_4_9(right_degree_4[9],d,d,d);
    cipher_matrix_jiang right_degree_4_10(right_degree_4[10],d,d,d);
    cipher_matrix_jiang right_degree_4_11(right_degree_4[11],d,d,d);
    cipher_matrix_jiang right_degree_4_12(right_degree_4[12],d,d,d);
    cipher_matrix_jiang right_degree_4_13(right_degree_4[13],d,d,d);
    cipher_matrix_jiang right_degree_4_14(right_degree_4[14],d,d,d);
    cipher_matrix_jiang right_degree_4_15(right_degree_4[15],d,d,d);


    cout<<" 4444444444444444444444444444"<<endl;
    // A^8
    vector<Ciphertext> left_degree_8=block_512_LeftSquare(context,encoder,scale,evaluator,relin_keys,gal_keys,
        left_degree_4_0,left_degree_4_1,left_degree_4_2,left_degree_4_3,left_degree_4_4,left_degree_4_5,
        left_degree_4_6,left_degree_4_7,left_degree_4_8,left_degree_4_9,left_degree_4_10,left_degree_4_11,
        left_degree_4_12,left_degree_4_13,left_degree_4_14,left_degree_4_15,left_degree_4_0,left_degree_4_1,left_degree_4_2,left_degree_4_3,left_degree_4_4,left_degree_4_5,
        left_degree_4_6,left_degree_4_7,left_degree_4_8,left_degree_4_9,left_degree_4_10,left_degree_4_11,
        left_degree_4_12,left_degree_4_13,left_degree_4_14,left_degree_4_15);
    //A+...+A^8
    vector<Ciphertext> right_degree_8=block_512_LeftMulRight_AddRight(context,encoder,scale,evaluator,relin_keys,gal_keys,
        left_degree_4_0,left_degree_4_1,left_degree_4_2,left_degree_4_3,left_degree_4_4,left_degree_4_5,
        left_degree_4_6,left_degree_4_7,left_degree_4_8,left_degree_4_9,left_degree_4_10,left_degree_4_11,
        left_degree_4_12,left_degree_4_13,left_degree_4_14,left_degree_4_15,right_degree_4_0,right_degree_4_1,right_degree_4_2,right_degree_4_3,right_degree_4_4,right_degree_4_5,
        right_degree_4_6,right_degree_4_7,right_degree_4_8,right_degree_4_9,right_degree_4_10,right_degree_4_11,
        right_degree_4_12,right_degree_4_13,right_degree_4_14,right_degree_4_15);
    
    cipher_matrix_jiang left_degree_8_0(left_degree_8[0],d,d,d);
    cipher_matrix_jiang left_degree_8_1(left_degree_8[1],d,d,d);
    cipher_matrix_jiang left_degree_8_2(left_degree_8[2],d,d,d);
    cipher_matrix_jiang left_degree_8_3(left_degree_8[3],d,d,d);
    cipher_matrix_jiang left_degree_8_4(left_degree_8[4],d,d,d);
    cipher_matrix_jiang left_degree_8_5(left_degree_8[5],d,d,d);
    cipher_matrix_jiang left_degree_8_6(left_degree_8[6],d,d,d);
    cipher_matrix_jiang left_degree_8_7(left_degree_8[7],d,d,d);
    cipher_matrix_jiang left_degree_8_8(left_degree_8[8],d,d,d);
    cipher_matrix_jiang left_degree_8_9(left_degree_8[9],d,d,d);
    cipher_matrix_jiang left_degree_8_10(left_degree_8[10],d,d,d);
    cipher_matrix_jiang left_degree_8_11(left_degree_8[11],d,d,d);
    cipher_matrix_jiang left_degree_8_12(left_degree_8[12],d,d,d);
    cipher_matrix_jiang left_degree_8_13(left_degree_8[13],d,d,d);
    cipher_matrix_jiang left_degree_8_14(left_degree_8[14],d,d,d);
    cipher_matrix_jiang left_degree_8_15(left_degree_8[15],d,d,d);

    cipher_matrix_jiang right_degree_8_0(right_degree_8[0],d,d,d);
    cipher_matrix_jiang right_degree_8_1(right_degree_8[1],d,d,d);
    cipher_matrix_jiang right_degree_8_2(right_degree_8[2],d,d,d);
    cipher_matrix_jiang right_degree_8_3(right_degree_8[3],d,d,d);
    cipher_matrix_jiang right_degree_8_4(right_degree_8[4],d,d,d);
    cipher_matrix_jiang right_degree_8_5(right_degree_8[5],d,d,d);
    cipher_matrix_jiang right_degree_8_6(right_degree_8[6],d,d,d);
    cipher_matrix_jiang right_degree_8_7(right_degree_8[7],d,d,d);
    cipher_matrix_jiang right_degree_8_8(right_degree_8[8],d,d,d);
    cipher_matrix_jiang right_degree_8_9(right_degree_8[9],d,d,d);
    cipher_matrix_jiang right_degree_8_10(right_degree_8[10],d,d,d);
    cipher_matrix_jiang right_degree_8_11(right_degree_8[11],d,d,d);
    cipher_matrix_jiang right_degree_8_12(right_degree_8[12],d,d,d);
    cipher_matrix_jiang right_degree_8_13(right_degree_8[13],d,d,d);
    cipher_matrix_jiang right_degree_8_14(right_degree_8[14],d,d,d);
    cipher_matrix_jiang right_degree_8_15(right_degree_8[15],d,d,d);
    cout<<" 88888888888888888888888888888888 "<<endl;
    // A^16
    vector<Ciphertext> left_degree_16=block_512_LeftSquare(context,encoder,scale,evaluator,relin_keys,gal_keys,
        left_degree_8_0,left_degree_8_1,left_degree_8_2,left_degree_8_3,left_degree_8_4,left_degree_8_5,
        left_degree_8_6,left_degree_8_7,left_degree_8_8,left_degree_8_9,left_degree_8_10,left_degree_8_11,
        left_degree_8_12,left_degree_8_13,left_degree_8_14,left_degree_8_15,left_degree_8_0,left_degree_8_1,left_degree_8_2,left_degree_8_3,left_degree_8_4,left_degree_8_5,
        left_degree_8_6,left_degree_8_7,left_degree_8_8,left_degree_8_9,left_degree_8_10,left_degree_8_11,
        left_degree_8_12,left_degree_8_13,left_degree_8_14,left_degree_8_15);
    //A+...+A^16
    vector<Ciphertext> right_degree_16=block_512_LeftMulRight_AddRight(context,encoder,scale,evaluator,relin_keys,gal_keys,
        left_degree_8_0,left_degree_8_1,left_degree_8_2,left_degree_8_3,left_degree_8_4,left_degree_8_5,
        left_degree_8_6,left_degree_8_7,left_degree_8_8,left_degree_8_9,left_degree_8_10,left_degree_8_11,
        left_degree_8_12,left_degree_8_13,left_degree_8_14,left_degree_8_15,right_degree_8_0,right_degree_8_1,
        right_degree_8_2,right_degree_8_3,right_degree_8_4,right_degree_8_5,
        right_degree_8_6,right_degree_8_7,right_degree_8_8,right_degree_8_9,right_degree_8_10,right_degree_8_11,
        right_degree_8_12,right_degree_8_13,right_degree_8_14,right_degree_8_15);
    
    cipher_matrix_jiang left_degree_16_0(left_degree_16[0],d,d,d);
    cipher_matrix_jiang left_degree_16_1(left_degree_16[1],d,d,d);
    cipher_matrix_jiang left_degree_16_2(left_degree_16[2],d,d,d);
    cipher_matrix_jiang left_degree_16_3(left_degree_16[3],d,d,d);
    cipher_matrix_jiang left_degree_16_4(left_degree_16[4],d,d,d);
    cipher_matrix_jiang left_degree_16_5(left_degree_16[5],d,d,d);
    cipher_matrix_jiang left_degree_16_6(left_degree_16[6],d,d,d);
    cipher_matrix_jiang left_degree_16_7(left_degree_16[7],d,d,d);
    cipher_matrix_jiang left_degree_16_8(left_degree_16[8],d,d,d);
    cipher_matrix_jiang left_degree_16_9(left_degree_16[9],d,d,d);
    cipher_matrix_jiang left_degree_16_10(left_degree_16[10],d,d,d);
    cipher_matrix_jiang left_degree_16_11(left_degree_16[11],d,d,d);
    cipher_matrix_jiang left_degree_16_12(left_degree_16[12],d,d,d);
    cipher_matrix_jiang left_degree_16_13(left_degree_16[13],d,d,d);
    cipher_matrix_jiang left_degree_16_14(left_degree_16[14],d,d,d);
    cipher_matrix_jiang left_degree_16_15(left_degree_16[15],d,d,d);

    cipher_matrix_jiang right_degree_16_0(right_degree_16[0],d,d,d);
    cipher_matrix_jiang right_degree_16_1(right_degree_16[1],d,d,d);
    cipher_matrix_jiang right_degree_16_2(right_degree_16[2],d,d,d);
    cipher_matrix_jiang right_degree_16_3(right_degree_16[3],d,d,d);
    cipher_matrix_jiang right_degree_16_4(right_degree_16[4],d,d,d);
    cipher_matrix_jiang right_degree_16_5(right_degree_16[5],d,d,d);
    cipher_matrix_jiang right_degree_16_6(right_degree_16[6],d,d,d);
    cipher_matrix_jiang right_degree_16_7(right_degree_16[7],d,d,d);
    cipher_matrix_jiang right_degree_16_8(right_degree_16[8],d,d,d);
    cipher_matrix_jiang right_degree_16_9(right_degree_16[9],d,d,d);
    cipher_matrix_jiang right_degree_16_10(right_degree_16[10],d,d,d);
    cipher_matrix_jiang right_degree_16_11(right_degree_16[11],d,d,d);
    cipher_matrix_jiang right_degree_16_12(right_degree_16[12],d,d,d);
    cipher_matrix_jiang right_degree_16_13(right_degree_16[13],d,d,d);
    cipher_matrix_jiang right_degree_16_14(right_degree_16[14],d,d,d);
    cipher_matrix_jiang right_degree_16_15(right_degree_16[15],d,d,d);



    Matrix<double> computed0,computed1,computed2,computed3,
    computed4,computed5,computed6,computed7,
    computed8,computed9,computed10,computed11,
    computed12,computed13,computed14,computed15;
    right_degree_16_0.dec_matrix_cipher(computed0,encoder,decryptor);
    right_degree_16_1.dec_matrix_cipher(computed1,encoder,decryptor);
    right_degree_16_2.dec_matrix_cipher(computed2,encoder,decryptor);
    right_degree_16_3.dec_matrix_cipher(computed3,encoder,decryptor);
    right_degree_16_4.dec_matrix_cipher(computed4,encoder,decryptor);
    right_degree_16_5.dec_matrix_cipher(computed5,encoder,decryptor);
    right_degree_16_6.dec_matrix_cipher(computed6,encoder,decryptor);
    right_degree_16_7.dec_matrix_cipher(computed7,encoder,decryptor);
    right_degree_16_8.dec_matrix_cipher(computed8,encoder,decryptor);
    right_degree_16_9.dec_matrix_cipher(computed9,encoder,decryptor);
    right_degree_16_10.dec_matrix_cipher(computed10,encoder,decryptor);
    right_degree_16_11.dec_matrix_cipher(computed11,encoder,decryptor);
    right_degree_16_12.dec_matrix_cipher(computed12,encoder,decryptor);
    right_degree_16_13.dec_matrix_cipher(computed13,encoder,decryptor);
    right_degree_16_14.dec_matrix_cipher(computed14,encoder,decryptor);
    right_degree_16_15.dec_matrix_cipher(computed15,encoder,decryptor);
    auto end_time = chrono::high_resolution_clock::now();
    auto time_diff = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << "The timing cost is: " << time_diff.count() / 1e6 << " s" << endl;
    // cout << "the computed C = A^2+A: " << endl;
    //  store_matrix(n,m,computed_C);
    // computed.print();
    std::ofstream outFile0("matrix_result_512_0.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile0 << computed0.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile0 << std::endl;
    }

    outFile0.close();
    std::ofstream outFile1("matrix_result_512_1.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile1 << computed1.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile1 << std::endl;
    }

    outFile1.close();
    std::ofstream outFile2("matrix_result_512_2.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile2 << computed2.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile2 << std::endl;
    }

    outFile2.close();
    std::ofstream outFile3("matrix_result_512_3.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile3 << computed3.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile3 << std::endl;
    }

    outFile3.close();
    std::ofstream outFile4("matrix_result_512_4.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile4 << computed4.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile4 << std::endl;
    }

    outFile4.close();
    std::ofstream outFile5("matrix_result_512_5.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile5 << computed5.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile5 << std::endl;
    }

    outFile5.close();
    std::ofstream outFile6("matrix_result_512_6.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile6 << computed6.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile6 << std::endl;
    }

    outFile6.close();
    std::ofstream outFile7("matrix_result_512_7.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile7 << computed7.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile7 << std::endl;
    }

    outFile7.close();
    std::ofstream outFile8("matrix_result_512_8.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile8 << computed8.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile8 << std::endl;
    }

    outFile8.close();
    std::ofstream outFile9("matrix_result_512_9.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile9 << computed9.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile9 << std::endl;
    }

    outFile9.close();
    std::ofstream outFile10("matrix_result_512_10.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile10 << computed10.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile10 << std::endl;
    }

    outFile10.close();
    std::ofstream outFile11("matrix_result_512_11.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile11 << computed11.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile11 << std::endl;
    }

    outFile11.close();
    std::ofstream outFile12("matrix_result_512_12.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile12 << computed12.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile12 << std::endl;
    }

    outFile12.close();
    std::ofstream outFile13("matrix_result_512_13.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile13 << computed13.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile13 << std::endl;
    }

    outFile13.close();
    std::ofstream outFile14("matrix_result_512_14.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile14 << computed14.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile14 << std::endl;
    }

    outFile14.close();
    std::ofstream outFile15("matrix_result_512_15.txt");
    for (size_t i = 0; i < d; i++) {
        for (size_t j = 0; j < d; j++) {
            outFile15 << computed15.get(j, i) << (j < n - 1 ? " " : "");
        }
        outFile15 << std::endl;
    }

    outFile15.close();


    cout<<"endl"<<"---------------------------------------------------------"<<endl;
    // cout << "the true C =  A^2+A: " << endl;
    // F.print();

}