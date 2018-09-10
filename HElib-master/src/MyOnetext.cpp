/* Copyright (C) 2012-2017 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */

/* Test_General.cpp - A general test program that uses a mix of operations over four ciphertexts.
 */
#include <NTL/ZZ.h>
#include <NTL/BasicThreadPool.h>
#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"
#include <NTL/lzz_pXFactoring.h>

#include <cassert>
#include <cstdio>

#ifdef DEBUG_PRINTOUT
#define debugCompare(ea,sk,p,c) {\
  NewPlaintextArray pp(ea);\
  ea.decrypt(c, sk, pp);\
  if (!equals(ea, pp, p)) { \
    std::cout << "oops:\n"; std::cout << p << "\n"; \
    std::cout << pp << "\n"; \
    exit(0); \
  }}
#else
#define debugCompare(ea,sk,p,c)
#endif

/**************

1. c1.multiplyBy(c0)
2. c0 += random constant
3. c2 *= random constant
4. tmp = c1
5. ea.shift(tmp, random amount in [-nSlots/2, nSlots/2])
6. c2 += tmp
7. ea.rotate(c2, random amount in [1-nSlots, nSlots-1])
8. c1.negate()
9. c3.multiplyBy(c2)
10. c0 -= c3

**************/
int main(int argc, char *argv[])
{
long m,p,r,L,c,w,d;
ZZX G;
r=2;
p=3;
m=0;
d=1;
L=1+NextPowerOfTwo(d);
if (m<2)
m = FindM(/*secprm=*/80, L, /*c=*/3, p, 1, 0, m, true);
FHEcontext context(m,p,r);


buildModChain(context,p,r);
cout<<context<<endl;
fstream keyFile("test.txt", fstream::out|fstream::trunc);
assert(keyFile.is_open());
//fstream keyFile1("test1.txt", fstream::out|fstream::trunc);
//assert(keyFile1.is_open());
G = context.alMod.getFactorsOverZZ()[0];
FHESecKey secretKey(context);
const FHEPubKey& publicKey=secretKey;
secretKey.GenSecKey(64);
addSome1DMatrices(secretKey);
EncryptedArray ea(context,G);
long nslots=ea.size();
PlaintextArray p0(ea);
p0.random();
p0.print(keyFile);
//cout<<p0<<endl;
Ctxt c0(publicKey);
ea.encrypt(c0,publicKey,p0);
//c0.print(keyFile1);
keyFile1<<c0;
PlaintextArray pp0(ea);
//ea.decrypt(c0,secretKey,pp0);
if(!pp0.equals(p0))
cerr<<"oops 0"<<endl;
else
cout<<"equal!"<<endl;
return 0;
}

