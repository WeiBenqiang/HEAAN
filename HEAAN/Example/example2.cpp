#include "../src/HEAAN.h"

using namespace std;
using namespace NTL;

int main()
{
  /*
  * Basic Parameters are in src/Params.h
  * If you want to use another parameter, you need cleto change src/Params.h file and re-complie this library.
  */
  // Parameters //
  long logq = 300; ///< Ciphertext modulus (this value should be <= logQ in "scr/Params.h")
  long logp = 30;  ///< Scaling Factor (larger logp will give you more accurate value)
  long logn = 3;   ///< number of slot is 1024 (this value should be < logN in "src/Params.h")
  long n = 1 << logn;
  long slots = n;
  long numThread = 1;

  // Construct and Generate Public Keys //
  cout << "Construct and Generate Public Keys" << endl;
  srand(time(NULL));
  //SetNumThreads(numThread);
  TimeUtils timeutils;
  Ring ring;
  SecretKey secretKey(ring);
  Scheme scheme(secretKey, ring);
  SchemeAlgo algo(scheme);

  // Make Random Array of Complex //
  //complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots);
  //complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(slots);

  // Make Random Array of real //
  //double *mvec1 = EvaluatorUtils::randomRealArray(slots);
  //double *mvec2 = EvaluatorUtils::randomRealArray(slots);

  //自己设置slots值
  double *mvec1 = new double[slots]{0.01, 0.2, 0.3, 0.4, 1.77, 1.85, 1.95, 1.99};
  double *mvec2 = new double[slots]{1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8};
  // Encrypt Two Arry of Complex //
  Ciphertext cipher1;
  scheme.encrypt(cipher1, mvec1, n, logp, logq);
  Ciphertext cipher2;
  scheme.encrypt(cipher2, mvec2, n, logp, logq);

  // double single = 0.2;
  Ciphertext ciphersingle, cinv, cpow, csig;
  // scheme.encryptSingle(ciphersingle, single, logp, logq);

  //function
  complex<double> *minv = new complex<double>[n];
  for (long i = 0; i < n; ++i)
  {
    minv[i] = 1. / mvec1[i];
  }

  //scheme.imultAndEqual(ciphersingle);
  //algo.inverse(cinv, cipher1, logp, 6); //(0<x<1)

  scheme.negateAndEqual(cipher1);
  scheme.addConstAndEqual(cipher1, 2, logp);

  complex<double> *decinv = scheme.decrypt(secretKey, cipher1);
  for (int i = 0; i < slots; ++i)
  {
    cout << "mvec1:: " << mvec1[i] << "   " << real(decinv[i]) << endl;
  }

  // Multiplication And Rescale //
  // Ciphertext cipherMult;
  // scheme.mult(cipherMult, cipher1, cipher2);
  // Ciphertext cipherMultAfterReScale;
  // scheme.reScaleBy(cipherMultAfterReScale, cipherMult, logp);

  //multByConst
  //Ciphertext cmultByConst;
  //scheme.multByConst(cmultByConst,cipher1,(complex<double>)5, 4);

  // Decrypt //
  //complex<double>* dvec1 = scheme.decrypt(secretKey, ciphersingle);
  //complex<double>* dvec2 = scheme.decrypt(secretKey, cipher2);

  //for(int i =0; i < n; i++)
  //    cout<<i<<" "<<mvec1[i]<<" "<<real(dvec1[i])<<endl;
  delete[] mvec1;
  delete[] mvec2;

  return 0;
}