#include "../src/HEAAN.h"
#include "../src/TimeUtils.h"
using namespace std;
using namespace NTL;

//求逆算法
void Inverse(Scheme &scheme, Ciphertext &cipher, long logp, long iters)
{
  Ciphertext ciphera, cipherb, cipherbPlus1;
  //初始值
  scheme.negate(ciphera, cipher);
  scheme.addConstAndEqual(ciphera, 2, logp);
  scheme.negate(cipherb, cipher);
  scheme.addConstAndEqual(cipherb, 1, logp);

  //迭代算法
  for (int i = 0; i < iters; ++i)
  {
    //b[n+1] = b[n]^2
    cout << "起始模数: " << ciphera.logq << " " << cipherb.logq << endl;
    scheme.squareAndEqual(cipherb);
    scheme.reScaleByAndEqual(cipherb, logp);

    //a[n+1] = a[n]*(1+b[n+1])
    scheme.addConst(cipherbPlus1, cipherb, 1, logp);
    cout << "before: " << ciphera.logq << " " << cipherbPlus1.logq << endl;
    long bitsDown = ciphera.logq - cipherbPlus1.logq;
    scheme.modDownByAndEqual(ciphera, bitsDown);
    cout << "after: " << ciphera.logq << " " << cipherbPlus1.logq << endl;
    scheme.multAndEqual(ciphera, cipherbPlus1);
    scheme.reScaleByAndEqual(ciphera, logp);

    cout << "结束模数: " << ciphera.logq << " " << cipherbPlus1.logq << endl;
  }
  cipher.copy(ciphera);
}

//平方根算法
void Sqrt(Scheme &scheme, Ciphertext &cipher, long logp, long iters)
{
  Ciphertext ciphera, cipherb, cipherbn;
  ciphera.copy(cipher);
  scheme.addConst(cipherb, cipher, -1, logp);

  //迭代算法、
  for (int i = 0; i < iters; ++i)
  {
    //第一个等式
    //a[n+1] = a[n]*(1 - 0.5 * b[n]
    // scheme.multByConst(cipherbn, cipherb, 0.5, logp);
    // scheme.reScaleByAndEqual(cipherbn, logp);

    scheme.divByPo2(cipherbn, cipherb, 1);
    scheme.negateAndEqual(cipherbn);
    scheme.addConstAndEqual(cipherbn, 1, logp);

    cout << "查看起始模数大小： a:" << ciphera.logq << " bn: " << cipherbn.logq << endl;
    long bitsDown = ciphera.logq - cipherbn.logq;
    scheme.modDownByAndEqual(ciphera, bitsDown);
    cout << "乘法之前模数大小: a" << ciphera.logq << " bn:" << cipherbn.logq << endl;

    scheme.multAndEqual(ciphera, cipherbn);
    scheme.reScaleByAndEqual(ciphera, logp);
    cout << "第一步之后模数大小： a" << ciphera.logq << "  bn:" << cipherbn.logq << endl;

    //第二个等式
    //b[n+1] = b[b]^2 * (b[n] - 3) / 4  更新cipherb
    Ciphertext cipherbSquare;
    scheme.square(cipherbSquare, cipherb);
    scheme.reScaleByAndEqual(cipherbSquare, logp); //760

    scheme.addConstAndEqual(cipherb, -3, logp);
    scheme.divByPo2AndEqual(cipherb, 2); //798
    cout << "第2步小括号之后： b_square: " << cipherbSquare.logq << "  b:" << cipherb.logq << endl;

    bitsDown = cipherb.logq - cipherbSquare.logq;
    scheme.modDownByAndEqual(cipherb, bitsDown);
    cout << "第2步乘法之前： b_square: " << cipherbSquare.logq << "  b:" << cipherb.logq << endl;

    scheme.multAndEqual(cipherb, cipherbSquare);
    scheme.reScaleByAndEqual(cipherb, logp);

    cout << "本轮结束之后： a" << ciphera.logq << "  b:" << cipherb.logq << endl;
    cout << "=========================================" << endl;
    //第一轮 b1= (mvec[i]-1)*(mvec[i]-1) * (mvec[i]-1 -3)/4
  }
  cipher.copy(ciphera);
}
//测试求逆算法
void testInverseByIterative(long logq, long logp, long logn, long iters)
{
  cout << "!!! START TEST InverseByIterative !!!" << endl;
  srand(time(NULL));
  SetNumThreads(8);
  TimeUtils timeutils;
  Ring ring;
  SecretKey secretKey(ring);
  Scheme scheme(secretKey, ring);

  long slots = (1 << logn);
  double *mvec = EvaluatorUtils::randomRealArray(slots);

  Ciphertext cipher;
  scheme.encrypt(cipher, mvec, slots, logp, logq);
  cout << "cipher logq before: " << cipher.logq << endl;

  timeutils.start("inversebyiterative");
  Inverse(scheme, cipher, logp, iters);
  timeutils.stop("inversebyiterative");

  complex<double> *dvec = scheme.decrypt(secretKey, cipher);
  for (int i = 0; i < slots; ++i)
  {
    cout << i << ": " << mvec[i] << " " << 1. / mvec[i] << " " << dvec[i] << endl;
  }
}

//测试求根号算法
void testSqrt(long logq, long logp, long logn, long iters)
{
  cout << "!!! START TEST Sqrt!!!" << endl;
  srand(time(NULL));
  SetNumThreads(8);
  TimeUtils timeutils;
  Ring ring;
  SecretKey secretKey(ring);
  Scheme scheme(secretKey, ring);

  long slots = (1 << logn);
  double *mvec = EvaluatorUtils::randomRealArray(slots);

  Ciphertext cipher, ciphera, cipherb, cipherbn;
  scheme.encrypt(cipher, mvec, slots, logp, logq);

  cout << "起始模数: " << ciphera.logq << " " << cipherb.logq << endl;

  timeutils.start("Sqrt");
  Sqrt(scheme, cipher, logp, iters);
  timeutils.stop("Sqrt");

  complex<double> *dvec = scheme.decrypt(secretKey, cipher);
  for (int i = 0; i < slots; ++i)
  {
    // cout << i << ": " << mvec[i] << " " <<(mvec[i]-1)*(mvec[i]-1) * (mvec[i]-1 -3)/4<< " " << dvec[i] << endl;
    cout << i << ": " << mvec[i] << " " << sqrt(mvec[i]) << " " << dvec[i] << endl;
  }
}

//测试求最大值算法
void testMAX(long logq, long logp, long logn, long iters)
{
  cout << "!!! START TEST MAX !!" << endl;
  srand(time(NULL));
  SetNumThreads(8);
  TimeUtils timeutils;
  Ring ring;
  SecretKey secretKey(ring);
  Scheme scheme(secretKey, ring);

  long slots = (1 << logn);
  double *mvec1 = EvaluatorUtils::randomRealArray(slots);
  double *mvec2 = EvaluatorUtils::randomRealArray(slots);
  double *max = new double[slots];
  //明文状态

  for (int i = 0; i < slots; ++i)
  {
    max[i] = (mvec1[i] > mvec2[i]) ? mvec1[i] : mvec2[i];
  }

  // x = (a+b)/2 ; y = (a-b) / 2
  // z = Sqrt(y^2)
  //Max(a,b) = x + z; Min(a,b) = x - z

  Ciphertext cipher1, cipher2, cipherX, cipherY;

  scheme.encrypt(cipher1, mvec1, slots, logp, logq);
  scheme.encrypt(cipher2, mvec2, slots, logp, logq);

  timeutils.start("MAX");
  scheme.add(cipherX, cipher1, cipher2);
  scheme.sub(cipherY, cipher1, cipher2);

  scheme.divByPo2AndEqual(cipherX, 1);
  scheme.divByPo2AndEqual(cipherY, 1);

  scheme.squareAndEqual(cipherY);
  scheme.reScaleByAndEqual(cipherY, logp);

  cout << "当前模数：cipherX: " << cipherX.logq << " cipherY: " << cipherY.logq << endl;

  Sqrt(scheme, cipherY, logp, iters);

  cout << "当前模数：cipherX: " << cipherX.logq << " cipherY: " << cipherY.logq << endl;

  long bitsDown = cipherX.logq - cipherY.logq;
  scheme.modDownByAndEqual(cipherX, bitsDown);

  scheme.addAndEqual(cipherX, cipherY);
  timeutils.stop("MAX");

  complex<double> *appMax = scheme.decrypt(secretKey, cipherX);

  // for (int i = 0; i < slots; ++i)
  // { //看看绝对值
  //   cout << i << ":" << mvec1[i] << " " << mvec2[i] << " true : " << abs((mvec1[i] - mvec2[i]) / 2)
  //        << "   not true:  " << dvec1[i] << endl;
  // }

  // for (int i = 0; i < slots; ++i)
  // { //看看最大值
  //   cout << i << ": " << mvec1[i] << " " << mvec2[i]
  //        << " true max: " << max[i] << "   approximate max:  " << appMax[i] << endl;
  // }

  delete[] max;
}

//======================================================================================
//=========        !!!!!!!!!!!!testCompare!!!!!!!!!!!!!!        ========================
//======================================================================================
void testCompare(long logq, long logp, long logn, long itersofDprime, long itersofD, long iters, long m)
{
  cout << "!!! START TEST COMPARE !!" << endl;
  srand(time(NULL));
  SetNumThreads(8);
  TimeUtils timeutils;
  Ring ring;
  SecretKey secretKey(ring);
  Scheme scheme(secretKey, ring);
  SchemeAlgo algo(scheme);

  long slots = (1 << logn);
  double *mvec1 = EvaluatorUtils::randomRealArray(slots);
  double *mvec2 = EvaluatorUtils::randomRealArray(slots);
  double *CompResult = new double[slots];

  //明文状态
  for (int i = 0; i < slots; ++i)
  {
    mvec1[i] = mvec1[i] + 0.5;
    mvec2[i] = mvec2[i] + 0.5;
    CompResult[i] = (mvec1[i] >= mvec2[i]) ? 1 : 0;
  }

  //加密
  Ciphertext cipher1, cipher2, cipherHalf, cipherPlus, ciphera, cipherb;

  scheme.encrypt(cipher1, mvec1, slots, logp, logq);
  scheme.encrypt(cipher2, mvec2, slots, logp, logq);

  timeutils.start("Compare");
  //------------------------------------------------------------------------
  //第一阶段  Normalize
  //a0 = (a/2) * Inv((a+b)./2 ;d')
  //------------------------------------------------------------------------
  scheme.add(cipherPlus, cipher1, cipher2);
  scheme.divByPo2AndEqual(cipherPlus, 1);

  Inverse(scheme, cipherPlus, logp, itersofDprime);
  cout << "inverse of cipherPlus modulus: " << cipherPlus.logq << endl;

  scheme.divByPo2(cipherHalf, cipher1, 1);
  cout << "inverse of cipherHalf modulus: " << cipherHalf.logq << endl;

  long bitsDown = cipherHalf.logq - cipherPlus.logq;
  scheme.modDownByAndEqual(cipherHalf, bitsDown);

  //乘积 求a0
  scheme.mult(ciphera, cipherHalf, cipherPlus);
  scheme.reScaleByAndEqual(ciphera, logp);
  //求b0
  scheme.negate(cipherb, ciphera);
  scheme.addConstAndEqual(cipherb, 1, logp);

  //------------------------------------------------------------------------==
  // 第二阶段
  //------------------------------------------------------------------------

  cout << "==============================================" << endl;
  Ciphertext cipheraPowerofm, cipherbPowerofm, cipherPlusofPower, cipherInvofPlus;
  for (int i = 0; i < iters; ++i) //m^iters = k
  {
    algo.power(cipheraPowerofm, ciphera, logp, m);
    algo.power(cipherbPowerofm, cipherb, logp, m);

    scheme.add(cipherPlusofPower, cipheraPowerofm, cipherbPowerofm);

    cout << "after power: cipheraPowerofm: " << cipheraPowerofm.logq
         << " cipheraPowerofm : " << cipherbPowerofm.logq << endl;
    //inv = Invese(a^m+b^m;d)
    cipherInvofPlus.copy(cipherPlusofPower);
    Inverse(scheme, cipherInvofPlus, logp, itersofD);

    cout << "求逆之后: cipheraPowerofm: " << cipheraPowerofm.logq
         << " cipherInvofPlus : " << cipherInvofPlus.logq << endl;

    //a[n+1] = a[n]^m * inv
    //b[n+1] = 1 -a[n+1]
    bitsDown = cipheraPowerofm.logq - cipherInvofPlus.logq;
    scheme.modDownByAndEqual(cipheraPowerofm, bitsDown);

    cout << "模切换之后: cipheraPowerofm: " << cipheraPowerofm.logq
         << " cipherInvofPlus : " << cipherInvofPlus.logq << endl;

    scheme.mult(ciphera, cipheraPowerofm, cipherInvofPlus);
    scheme.reScaleByAndEqual(ciphera, logp);

    scheme.negate(cipherb, ciphera);
    scheme.addConstAndEqual(cipherb, 1, logp);

    cout << "==============================================" << endl;
  }
  timeutils.stop("Compare");

  complex<double> *appComp = scheme.decrypt(secretKey, ciphera);

  for (int i = 0; i < slots; ++i)
  {
    cout << i << ": " << mvec1[i] << " " << mvec2[i]
         << " true value: " << CompResult[i] //<< endl;
         << "   approximate max:  " << appComp[i] << endl;
  }

  delete[] CompResult;
}

//======================================================================================
//=========        !!!!!!!!!!!!testEfficientCompare!!!!!!!!!!!!!!        ===============
//======================================================================================

void testEfficientCompare(long logq, long logp, long logn, int degreeOff, long itersOff)
{
  cout << "!!! START TEST EfficientCompare !!" << endl;
  srand(time(NULL));
  SetNumThreads(8);
  TimeUtils timeutils;
  Ring ring;
  SecretKey secretKey(ring);
  Scheme scheme(secretKey, ring);
  SchemeAlgo algo(scheme);

  long slots = (1 << logn);
  double *mvec1 = EvaluatorUtils::randomRealArray(slots);
  double *mvec2 = EvaluatorUtils::randomRealArray(slots);
  double *CompResult = new double[slots];

  //明文状态
  for (int i = 0; i < slots; ++i)
  {
    mvec1[i] = mvec1[i] + 0.5;
    mvec2[i] = mvec2[i] + 0.5;
    //CompResult[i] = mvec1[i] - mvec2[i];
    CompResult[i] = (mvec1[i] >= mvec2[i]) ? 1 : 0;
  }

  //加密
  Ciphertext cipher1, cipher2, cipherX, res;

  scheme.encrypt(cipher1, mvec1, slots, logp, logq);
  scheme.encrypt(cipher2, mvec2, slots, logp, logq);

  timeutils.start("EfficientCompare");
  //x = a-b;
  scheme.sub(cipherX, cipher1, cipher2);

  double *coeffsOff;
  if (degreeOff == 3)
  {
    coeffsOff = new double[4]{0, 3. / 2, 0, -1. / 2};
  }
  else if (degreeOff == 5)
  {
    coeffsOff = new double[6]{0, 15. / 8, 0, -10. / 8, 0, 3. / 8};
  }
  else if (degreeOff == 7)
  {
    coeffsOff = new double[8]{0, 35. / 16, 0, -35. / 16, 0, 21. / 16, 0, -5. / 16};
  }
  else if (degreeOff == 9)
  {
    coeffsOff = new double[10]{0, 315. / 128, 0, -420. / 128, 0, 378. / 128,
                               0, -180. / 128, 0, 35. / 128};
  }

  long dlogp = 2 * logp;
  Ciphertext *cpows = new Ciphertext[degreeOff];

  for (int i = 0; i < itersOff; ++i)
  {
    // 构造 x 的 pows 数组
    algo.powerExtended(cpows, cipherX, logp, degreeOff);
    //cout << "powof " << i + 1 << "的模数 " << cpows[i].logq << endl;

    scheme.multByConst(res, cpows[0], coeffsOff[1], logp);
    scheme.addConstAndEqual(res, coeffsOff[0], dlogp);

    Ciphertext aixi;
    for (int i = 1; i < degreeOff; ++i)
    {
      if (abs(coeffsOff[i + 1]) > 1e-27)
      {
        scheme.multByConst(aixi, cpows[i], coeffsOff[i + 1], logp);
        scheme.modDownToAndEqual(res, aixi.logq);
        scheme.addAndEqual(res, aixi);
      }
    }
    scheme.reScaleByAndEqual(res, logp);
    cipherX.copy(res);

    cout << "iters = " << i << " ciherX的模数：" << cipherX.logq << endl;
    cout << "==============================================" << endl;
  }

  // return (x+1)/2
  scheme.addConstAndEqual(cipherX, 1, logp);
  scheme.divByPo2AndEqual(cipherX, 1);
  timeutils.stop("EfficientCompare");

  complex<double> *appComp = scheme.decrypt(secretKey, cipherX);

  for (int i = 0; i < slots; ++i)
  {
    cout << i << ": " << mvec1[i] << " " << mvec2[i]
         << " true value: " << CompResult[i]
         //<< -0.5 * pow(CompResult[i], 3) + 1.5 * CompResult[i]
         //  << (3. / 8 * pow(CompResult[i], 5)
         //  - 10. / 8 * pow(CompResult[i], 3)
         //  + 15. / 8 * CompResult[i] + 1) /2
         << "   approximate max:  " << appComp[i] << endl;
  }

  delete[] CompResult;
}

//======================================================================================
//=========        !!!!!!!!!!!!testMoreEfficientCompare!!!!!!!!!!!!!!        ===========
//======================================================================================
void testMoreEfficientCompare(long logq, long logp, long logn,
                              int degreeOff, int degreeOfg, long itersOff, long itersOfg)
{
  cout << "!!! START TEST EfficientCompare !!" << endl;
  srand(time(NULL));
  //SetNumThreads(8);
  TimeUtils timeutils;
  Ring ring;
  SecretKey secretKey(ring);
  Scheme scheme(secretKey, ring);
  SchemeAlgo algo(scheme);

  long slots = (1 << logn);
  double *mvec1 = EvaluatorUtils::randomRealArray(slots);
  double *mvec2 = EvaluatorUtils::randomRealArray(slots);
  double *CompResult = new double[slots];

  //明文状态
  for (int i = 0; i < slots; ++i)
  {
    mvec1[i] = mvec1[i] + 0.5;
    mvec2[i] = mvec2[i] + 0.5;
    //CompResult[i] = mvec1[i] - mvec2[i];
    CompResult[i] = (mvec1[i] >= mvec2[i]) ? 1 : 0;
  }

  //加密
  Ciphertext cipher1, cipher2, cipherX, res;

  scheme.encrypt(cipher1, mvec1, slots, logp, logq);
  scheme.encrypt(cipher2, mvec2, slots, logp, logq);

  timeutils.start("MoreEfficientCompare");
  //x = a-b;
  scheme.sub(cipherX, cipher1, cipher2);

  double *coeffsOff, *coeffsOfg;

  // g的循环 (t = 1/4)
  if (degreeOfg == 3)
  {
    coeffsOfg = new double[4]{0, 2126. / 1024, 0, -1359. / 1024};
  }
  else if (degreeOfg == 5)
  {
    coeffsOfg = new double[6]{0, 3334. / 1024, 0, -6108. / 1024, 0, 3796. / 1024};
  }
  else if (degreeOfg == 7)
  {
    coeffsOfg = new double[8]{0, 4589. / 1024, 0, -16577. / 1024, 0,
                              25614. / 1024, 0, -12860. / 1024};
  }
  else if (degreeOfg == 9)
  {
    coeffsOfg = new double[10]{0, 5850. / 1024, 0, -34974. / 1024, 0, 97015. / 1024,
                               0, -113492. / 1024, 0, 46623. / 1024};
  }
  long dlogp = 2 * logp;
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // g的循环
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
  Ciphertext *cpowsOfg = new Ciphertext[degreeOfg];

  for (int i = 0; i < itersOfg; ++i)
  {
    // 构造 x 的 pows 数组
    algo.powerExtended(cpowsOfg, cipherX, logp, degreeOfg);
    //cout << "powof " << i + 1 << "的模数 " << cpows[i].logq << endl;

    scheme.multByConst(res, cpowsOfg[0], coeffsOfg[1], logp);
    scheme.addConstAndEqual(res, coeffsOfg[0], dlogp);

    Ciphertext aixi;
    for (int i = 1; i < degreeOfg; ++i)
    {
      if (abs(coeffsOfg[i + 1]) > 1e-27)
      {
        scheme.multByConst(aixi, cpowsOfg[i], coeffsOfg[i + 1], logp);
        scheme.modDownToAndEqual(res, aixi.logq);
        scheme.addAndEqual(res, aixi);
      }
    }
    scheme.reScaleByAndEqual(res, logp);
    cipherX.copy(res);

    cout << "itersOf f = " << i << " ciherX的模数：" << cipherX.logq << endl;
    cout << "==============================================" << endl;
  }

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // f的循环
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
  if (degreeOff == 3)
  {
    coeffsOff = new double[4]{0, 3. / 2, 0, -1. / 2};
  }
  else if (degreeOff == 5)
  {
    coeffsOff = new double[6]{0, 15. / 8, 0, -10. / 8, 0, 3. / 8};
  }
  else if (degreeOff == 7)
  {
    coeffsOff = new double[8]{0, 35. / 16, 0, -35. / 16, 0, 21. / 16, 0, -5. / 16};
  }
  else if (degreeOff == 9)
  {
    coeffsOff = new double[10]{0, 315. / 128, 0, -420. / 128, 0, 378. / 128,
                               0, -180. / 128, 0, 35. / 128};
  }

  Ciphertext *cpowsOff = new Ciphertext[degreeOff];

  for (int i = 0; i < itersOff; ++i)
  {
    // 构造 x 的 pows 数组
    algo.powerExtended(cpowsOff, cipherX, logp, degreeOff);
    //cout << "powof " << i + 1 << "的模数 " << cpows[i].logq << endl;

    scheme.multByConst(res, cpowsOff[0], coeffsOff[1], logp);
    scheme.addConstAndEqual(res, coeffsOff[0], dlogp);

    Ciphertext aixi;
    for (int i = 1; i < degreeOff; ++i)
    {
      if (abs(coeffsOff[i + 1]) > 1e-27)
      {
        scheme.multByConst(aixi, cpowsOff[i], coeffsOff[i + 1], logp);
        scheme.modDownToAndEqual(res, aixi.logq);
        scheme.addAndEqual(res, aixi);
      }
    }
    scheme.reScaleByAndEqual(res, logp);
    cipherX.copy(res);

    cout << "itersOf g = " << i << " ciherX的模数：" << cipherX.logq << endl;
    cout << "==============================================" << endl;
  }

  //-------------------------------------------------------------------
  // return (x+1)/2
  //-------------------------------------------------------------------
  scheme.addConstAndEqual(cipherX, 1, logp);
  scheme.divByPo2AndEqual(cipherX, 1);

  timeutils.stop("MoreEfficientCompare");

  complex<double> *appComp = scheme.decrypt(secretKey, cipherX);

  for (int i = 0; i < slots; ++i)
  {
    cout << i << ": " << mvec1[i] << " " << mvec2[i]
         << " true value: " << CompResult[i]
         << "   approximate max:  " << appComp[i] << endl;
  }

  delete[] CompResult;
}

//============================================================================================
int main()
{
  /*
  * Basic Parameters are in src/Params.h
  * If you want to use another parameter, you need to change src/Params.h file and re-complie this library.
  */

  // Parameters //
  long logq = 800; ///< Ciphertext modulus (this value should be <= logQ in "scr/Params.h")
  long logp = 40;  ///< Scaling Factor (larger logp will give you more accurate value)
  long logn = 10;   ///< number of slot is 1024 (this value should be < logN in "src/Params.h")

  long itersOfInverse = 5;
  long itersOfSqrt = 6;
  long itersOfMax = 11;
  long itersOfCompDprime = 4;
  long itersOfCompD = 4;
  long itersOfComp = 1;
  long m = 3; // m^t = k  2^3 = 8

  long itersofEfficientComp = 3;

  int degreeOff = 3;

  int degreeOfg = 3;
  long itersOff = 3;
  long itersOfg = 2;

  //long numThread = 1;

  //my test
  //testInverseByIterative(logq, logp, logn, itersOfInverse);
  //testSqrt(logq, logp, logn, itersOfSqrt);

  testMAX(logq, logp, logn, itersOfMax); //求最大值
  //testCompare(logq, logp, logn, itersOfCompDprime, itersOfCompD, itersOfComp, m); //求比较

  // testEfficientCompare(logq, logp, logn, degreeOff, itersofEfficientComp);

  //testMoreEfficientCompare(logq, logp, logn,
  //                         degreeOff, degreeOfg, itersOff, itersOfg);
  return 0;
}