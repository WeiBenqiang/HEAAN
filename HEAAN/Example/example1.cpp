#include "../src/HEAAN.h"

using namespace std;
using namespace NTL;

int main() {
  /*
  * Basic Parameters are in src/Params.h
  * If you want to use another parameter, you need cleto change src/Params.h file and re-complie this library.
  */

  // Parameters //
  long logq = 300; ///< Ciphertext modulus (this value should be <= logQ in "scr/Params.h")
  long logp = 30; ///< Scaling Factor (larger logp will give you more accurate value)
  long logn = 3; ///< number of slot is 1024 (this value should be < logN in "src/Params.h")
  long n = 1 << logn;
  long slots = n;
  long numThread = 1;
	
  // Construct and Generate Public Keys //
  cout<<"Construct and Generate Public Keys"<<endl;
  srand(time(NULL));
  //SetNumThreads(numThread);
  TimeUtils timeutils;
  Ring ring;
  SecretKey secretKey(ring);
  Scheme scheme(secretKey, ring);

  // Make Random Array of Complex //
  //complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots);
  //complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(slots);
  
  // Make Random Array of real //
  double *mvec1 = EvaluatorUtils::randomRealArray(slots);
  double *mvec2 = EvaluatorUtils::randomRealArray(slots);

  // Encrypt Two Arry of Complex //
  Ciphertext cipher1;
  scheme.encrypt(cipher1, mvec1, n, logp, logq);
  Ciphertext cipher2;
  scheme.encrypt(cipher2, mvec2, n, logp, logq);
  
  // Addition //
  Ciphertext cipherAdd;
  timeutils.start("add");
  scheme.add(cipherAdd, cipher1, cipher2);
  timeutils.stop("add");

  // Multiplication And Rescale //
  Ciphertext cipherMult;
  scheme.mult(cipherMult, cipher1, cipher2);

  Ciphertext cipherMultAfterReScale;
  scheme.reScaleBy(cipherMultAfterReScale, cipherMult, logp);
  
  // Rotation //
  long idx = 1;
  Ciphertext cipherRot;
  scheme.addLeftRotKey(secretKey,idx); ///< When you need left rotation for the vectorized message
  scheme.leftRotateFast(cipherRot, cipher1, idx);
  
  // Decrypt //
  complex<double>* dvec1 = scheme.decrypt(secretKey, cipher1);

  complex<double>* dvec2 = scheme.decrypt(secretKey, cipher2);
  //StringUtils::compare(mvec1, (double *) real(dvec1), n, "encrypt");
  for(int i =0; i < n; i++)
      cout<<i<<" "<<mvec1[i]<<" "<<real(dvec1[i])<<endl;
  
  return 0;

}