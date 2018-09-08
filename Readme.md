##Transform_based_Tensor_Decomposition_FHE

###1、The core source code of the thesis is described as follows:"PaperCode.m"

Our experimental data is video. Aerial data (aerial data) and perception data sensory data (Intelligent building and Intelligent Transportation case)

All homomorphic encryption and decryption experiments are involved.

You need to set the HElib library in the Linux configuration and then test the code of the HElib.

Involved basic code S-tSVD algorithm

#### 1.1 Basic Code S-tSVD Algorithm

obj = VideoReader('Aerial_data.mp4');

numFrames = obj.NumberOfFrames; % Total number of frames 7495

n=3; %This can be set by yourself.

D=ones(110,128,n); %Create an empty 3D matrix.

i=1;

array=[];  %This array is defined to record the value of the captured frame.

for k = 500: 2500: numFrames % If the number of captured frames is the same as that of n, 500: 2500 1:49

​    array(i)=k; 

​    frame = read(obj,k);%Indicates the number of the frame to be read.

​     %Every time you read a picture, convert it into a gray picture.

​    A=rgb2gray(frame); 

​    J=imresize(A,[110,128]);

​    A1=double(J);

​    D(:,:,i)=A1;

​    i=i+1;​

% tSVD decomposition

x_hat=fft(D,[],3);

for i=1:size(D,n)

   [U1,S1,V1]=svd(x_hat(:,:,i));

​    if i==1

​       D1=ones(size(U1,1),size(U1,2),size(D,3));

​       D2=ones(size(S1,1),size(S1,2),size(D,3));

​       D3=ones(size(V1,1),size(V1,2),size(D,3));

​    end

   D1(:,:,i) = U1;

   D2(:,:,i) = S1;

   D3(:,:,i) = V1;

end

U=ifft(D1,[],3);

S=ifft(D2,[],3);

V=ifft(D3,[],3);   

#### 1.2 Calculate compression degree (DimensionalityReduction Ratio)

%---Estimation Approximation Degree (Dimension Reduction Experiment)

% Select the original __ of the first __ non-zero singular values

​    for n1 = 1:5 %（1）

​       for n2 = 1:5 %（2）

​           for n3=1:3

​           S_1(n1,n2,n3) = S(n1,n2,n3);  

​           end

​       end 

​    end 

​    %Left singular tensor

​    for n1 = 1:110 

​       for n2 = 1:5  %（3）

​           for n3 = 1:3

​           U_1(n1,n2,n3) = U(n1,n2,n3);  

​           end

​       end 

​    end  

​    %Right singular tensor

​    for n1 = 1:128

​       for n2 = 1:5 %（4）

​           for n3= 1:3

​           V_1(n1,n2,n3) = V(n1,n2,n3);  

​           end

​       end 

​    end  

   %Calculate the product  New Tensor Multiplication of the tensor.

%display (size (S_1));

c =t_product (U_1, S_1); %c is calculated as a tensor.

% V_1 transposition

​    h2= ones(5,128,3); %（5）-1Location

​    for j=1:size(D,n)

​       h2(:,:,j)=V_1(:,:,j)'; %110*126*3

​    end

​    for k=2:size(D,n)

​       t1=h2(:,:,k);

​       h2(:,:,k)=h2(:,:,2+2-k);

​       h2(:,:,2+2-k)=t1;

​    end

   c =t_product(c,h2);

  result=[c(:,:,1);c(:,:,2);c(:,:,3)];

   %imshow(uint8(result));

​    %Relative error

   T_error=result-B1;

   display(T_error);

​    %Calculate the error rate of reconstruction

​    n =norm(T_error, 'fro' );

​    %display(n);

​    %Frobenius norm of the original tensor

   n_B=norm(B1,'fro');

​    %Tensor ApproximationRatio(reconstruction error degree)

   e=n/n_B; 

   display(e);

   RES=20*log10(e);

   display(RES);

​    p =size(D,1)*size(D,2)*size(D,3)/(n2*(size(D,1)+size(D,2)+1));

   display(1/p);

#### 1.3 Comparison with traditional svd and tSVD-slice

​    %traditional svd method

   addpath /Users/apple/Downloads/tensor_toolbox 

​    obj= VideoReader('Aerial_data.mp4');

   numFrames = obj.NumberOfFrames;% 7495

   D=ones(110,128,3); %Create an empty 3D matrix

   i=1;

​    for k = 500 : 2500 : numFrames

​        frame = read(obj,k);%Indicates the number of the frame to be read.

Every time you read a picture, convert him into a gray picture

​        A=rgb2gray(frame); 

​        J=imresize(A,[110,128]);

​        A1=double(J);

​        D(:,:,i)=A1;

​        i=i+1;

​    end

   A=tensor(D);

​    tic

   A1=tenmat(A,1);  %Expand mode-1 to n1* (n2*n3).

%A1=A1 ';

[U, S, V]=svd (A1.data);

%toc

% Select the first __ non-zero singular values and multiply them

​    %tic

   U1=U(:,1:53); %(1)Modify

   S1=S(1:53,1:53); %(2)Modify

   V1=V(:,1:53); %(3)Modify

   AA=U1*S1'*V1';

​    toc

   A1=double(A1);  

   SVD_Error = AA-A1;

   n_svd = norm(SVD_Error, 'fro' );

   n_B_svd = norm(A1,'fro');

   e_svd = n_svd/n_B_svd; 

   display(e_svd);

   RES_svd=20*log10(e_svd);

   display(RES_svd);

​    %--------------------

​    %tSVD-tubal method

​    obj= VideoReader('Aerial_data.mp4');

   numFrames = obj.NumberOfFrames;

   n=3;

   D=ones(110,128,n);%Create an empty 3D matrix

   i=1;

​    for k = 500 : 2500 : numFrames

​        frame = read(obj,k);%Indicates the number of the frame to be read

​        A=rgb2gray(frame); 

​        J=imresize(A,[110,128]);

​        A1=double(J);

​        D(:,:,i)=A1;

​        i=i+1;

​    end

​    tic %t-product algorithm without fft frequency domain change

​     for i=1:size(D,n)

​       [U1,S1,V1]=svd(D(:,:,i));

​       if i==1

​           D1=ones(size(U1,1),size(U1,2),size(D,n));

​           D2=ones(size(S1,1),size(S1,2),size(D,n));

​           D3=ones(size(V1,1),size(V1,2),size(D,n));

​       end

​       D1(:,:,i) = U1;

​       D2(:,:,i) = S1;

​       D3(:,:,i) = V1;

​     end

​     %Implement dimension reduction

​     for n1 = 1:35 %（1）

​       for n2 = 1:35 %（2）

​           for n3=1:3

​           S_1(n1,n2,n3) = D2(n1,n2,n3);  

​           end

​       end 

​    end 

​    %Left singular tensor

​    for n1 = 1:110 

​       for n2 = 1:35 %（3）

​           for n3 = 1:3

​           U_1(n1,n2,n3) = D1(n1,n2,n3); 

​           end

​       end 

​    end  

​    %Right singular tensor

​    for n1 = 1:128

​       for n2 = 1:35 %（4）

​           for n3=1:3

​           V_1(n1,n2,n3) = D3(n1,n2,n3);

​           end

​       end 

​    end  

%  Matrix multiplication

​    for m = 1:1:3

​       mid(:,:,m) = U_1(:,:,m)*S_1(:,:,m);

​    end

​    for n=1:1:3

​       result1(:,:,n)= mid(:,:,n)*V_1(:,:,n)';

​    end

​    toc

​    %new tensor product

%    c = t_product(U_1,S_1); %cA tensor is calculated

%    % V_1Transposition

%    h2 = ones(35,128,3); %（5）-1Location

%    for j=1:size(D,n)

%         h2(:,:,j)=V_1(:,:,j)'; %110*126*3

%    end

%    for k=2:size(D,n)

%         t1=h2(:,:,k);

%         h2(:,:,k)=h2(:,:,2+2-k);

%         h2(:,:,2+2-k)=t1;

%    end

%   c = t_product(c,h2);

%   result=[c(:,:,1);c(:,:,2);c(:,:,3)];

​    %-------Result presentation----------

   B1=[D(:,:,2);D(:,:,3);D(:,:,1)];

​    %imshow(uint8(B1));

​    %imshow(uint8(D));

   result=[result1(:,:,2);result1(:,:,3);result1(:,:,1)];

​    %imshow(uint8(result));

​    %Relative error

   T_error=result-B1;

​    %display(T_error);

​    %Calculate the error rate of reconstruction

​    n =norm(T_error, 'fro' );

   n_B=norm(B1,'fro');

   e=n/n_B; 

   display(e);

   RES=20*log10(e);

   display(RES);

#### 1.4 Involved Computing Running Timerunning time comparetraditional svd

​    %traditional svd time

   addpath /Users/apple/Downloads/tensor_toolbox 

​    obj= VideoReader('Aerial_data.mp4');

   numFrames = obj.NumberOfFrames; 

​    D=ones(100,200,30);%Create an empty 3D matrix% first where you want to modify

   i=1;

​    for k = 500 : 2500 : numFrames%

​        frame = read(obj,k);%Indicates the number of the frame to be read.

Every time you read a picture, convert him into a gray picture

​        A=rgb2gray(frame); 

​        J=imresize(A,[100,200]); %The second place to be modified

​        A1=double(J);

​        D(:,:,i)=A1;

​        i=i+1;

​    end

   A=tensor(D);

​    tic

   A1=tenmat(A,1);  

   [U,S,V]=svd(A1.data); 

   toc  %[1,2,3]:0.016323  [10,20,3]:0.026745 [10,20,30]:0.027036  [10,200,30]:0.701844[100,200,30]:5.219157  [100,200,300]Errorusing svd

​    %Requested 60000x60000 (26.8GB)array exceeds maximum array size preference. Creation of arrays greater than

​    %this limit may take a longtime and cause MATLAB to become unresponsive. See array size limit orpreference

​    %panel for more information.

​    %S-tSVD

​    obj= VideoReader('Aerial_data.mp4');

   numFrames = obj.NumberOfFrames;% 9243 

   D=ones(10,200,30);%Create an empty 3D matrix

   i=1;

​    for k = 1 : 309: numFrames% 3081 309 31

​        frame = read(obj,k);%Indicates the number of the frame to be read.

Every time you read a picture, convert him into a gray picture

​        A=rgb2gray(frame); 

​        J=imresize(A,[10,200]);

​        A1=double(J);

​        D(:,:,i)=A1;

​        i=i+1;

​    end

   x_hat=fft(D,[],3); %Any cyclic matrix can be diagonalized by Fourier transform matrix

​    tic

​    for i=1:size(D,3)

​       [U1,S1,V1]=svd(x_hat(:,:,i));

​       if i==1

​           D1=ones(size(U1,1),size(U1,2),size(D,3));

​           D2=ones(size(S1,1),size(S1,2),size(D,3));

​           D3=ones(size(V1,1),size(V1,2),size(D,3));

​       end

​       D1(:,:,i) = U1;

​       D2(:,:,i) = S1;

​       D3(:,:,i) = V1;

​    end

​    toc%[1,2,3]:0.010292[10,20,3]:0.010729 [10,20,30]:0.015893 [100,200,30]:0.353426[100,200,300]:1.760767 [1000,2000,3000]:22.595927

   U=ifft(D1,[],3);

   S=ifft(D2,[],3);

   V=ifft(D3,[],3);

#### 1.5 The feature matrix and eigenvalues obtained after tensor decomposition are obtained

   max_S2=max(max(S(:,:,2)));

   display(max_S2);

   min_S22=min(min(S(:,:,2)));

   display(min_S22);

​    %S（：，；，k）Information about each slice

   a1=1:1:50;

   b1=1:1:50;

   [x1,y1]=meshgrid(a1,b1);

   M_S001=ones(50,50); %110*128

​    %display(S(1,1,1));

​    for j=1:50

​       for k=1:50

​       M_S001(k,j)=S(k,j,21);

​       end

​    end

   p001=mesh(x1,y1,M_S001);

   c=mean(M_S001(:)); %Average

   display(c)

   title('S(:,:,21)')

​    %S（：，k，：）Slice information

​    T_S2=S(:,2,:);

​    display(T_S2); %110*1*155

​    M_S2=ones(110,155);

​    %  for i=1:155

​    %       M_S2(:,i)=T_S2(:,:,i); %110行1列

​    %  end

​     for i=1:110

​        for j=1:155

​            M_S2(i,j)=T_S2(i,1,j);

​        end

​     end

​    c2=mean(M_S2(:)); %Average

​    display(c2); 

​    p2=mesh(x,y,M_S2);

​    title('S(:,2,:)');

​     %S（：，j，：）Slice information

​    S_matrix1=S(:,1,:);

​    S_matrix2=S(:,2,:);

​    S_matrix3=S(:,3,:);

​    S_matrix4=S(:,4,:);

​    S_vactor1=S_matrix1(:); %42240*1

​    S_vactor2=S_matrix2(:);

​    S_vactor3=S_matrix3(:);

​    S_vactor4=S_matrix4(:);

​    I_index_number1=find(S_vactor1~=0); %..*1

​    I_index_number2=find(S_vactor2~=0);

​    I_index_number3=find(S_vactor3~=0);

​    I_index_number4=find(S_vactor4~=0);

​    x=1:155;

​    plot(x,S_vactor1(I_index_number1,1),'.-r','MarkerSize',10);

​    hold on

​    plot(x,S_vactor2(I_index_number2,1),'.-b','MarkerSize',10);

​    hold on

​    plot(x,S_vactor3(I_index_number3,1),'.-k','MarkerSize',10);

​    hold on

​    plot(x,S_vactor4(I_index_number4,1),'.-g','MarkerSize',10);

​    xlabel('Experiment');

​    ylabel('SingularValue');

​     %set(gca,'xTick',1:160);

​    legend('Slice one','Slice two','Slice three','Slice four');

###2、new tensor product (t-product)  "t_production.m"

function c=t_product(a,b)

%     a_n1=110;a_n2=110;a_n3=3;

%     b_n1=110;b_n2=110;b_n3=3;

   [a_n1 a_n2 a_n3]=size(a); %n1 n2 n3 values 110 110 3

   [b_n1 b_n2 b_n3]=size(b); %110 110 3

   c=zeros(a_n1,b_n2,a_n3); 

   A=cell(a_n3,1);

   B=cell(b_n3,1);

​    for i=1:a_n3

​        A{i}=a(:,:,i);

​       B{i}=b(:,:,i);

​    end

   index_up=zeros(1,a_n3);

   index_down=zeros(1,a_n3);

​    for i=1:a_n3   

​       index_up(i)=a_n3-i+1;

​       index_down(i)=i;

​    end

​    %sStores the left items of the expanded tensor

   s=cell(a_n3,a_n3);

​    for i=1:a_n3

​       for j=1:a_n3

​           if i==j

​                s{i,j}=A{1};

​           end      

​           if j>i

​                s{i,j}=A{index_up(j-i)};

​           end      

​           if j<i

​                s{i,j}=A{index_down(i-j+1)};

​           end     

​       end  

​    end

​    %Apply for a memory size for storing calculation results

   re=cell(a_n3,1);

​    for i=1:a_n3

​       re{i}=zeros(a_n1,b_n2);

​    end

%    display(re);

​    for i=1:a_n3

​       for j=1:a_n3

​           for k=1:1

​                re{i,k}=re{i,k}+s{i,j}*B{j,k};

​           end

​       end   

​    end

% Calculate the solution of the product

​    for i=1:a_n3

​       c(:,:,i)=re{i};        

​    end

end

###3、Code for full homomorphic encryption based on the HElib library

/* Copyright (C) 2012-2017 IBMCorp.

 * This program is Licensed under the ApacheLicense, Version 2.0

 * (the "License"); you may not usethis file except in compliance

 * with the License. You may obtain a copy ofthe License at

 *  http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreedto in writing, software

 * distributed under the License is distributedon an "AS IS" BASIS,

 * WITHOUT WARRANTIES OR CONDITIONS OF ANYKIND, either express or implied.

 * See the License for the specific languagegoverning permissions and

 * limitations under the License. Seeaccompanying LICENSE file.

 */

/* Test_General.cpp - Ageneral test program that uses a mix of operations over four ciphertexts.

 */

\#include <NTL/ZZ.h>

\#include<NTL/BasicThreadPool.h>

\#include "FHE.h"

\#include "timing.h"

\#include"EncryptedArray.h"

\#include <NTL/lzz_pXFactoring.h>

\#include "Ctxt.h"

\#include <stdio.h>

\#include <cassert>

\#include <cstdio>

 

//#ifdef DEBUG_PRINTOUT

//#definedebugCompare(ea,sk,p,c) {\

//  NewPlaintextArray pp(ea);\

  ea.decrypt(c, sk, pp);\

  if (!equals(ea, pp, p)) { \

​    std::cout << "oops:\n";std::cout << p << "\n"; \

​    std::cout << pp <<"\n"; \

​    exit(0); \

  }}

//#else

//#definedebugCompare(ea,sk,p,c)

//#endif

 

/**************

\1. c1.multiplyBy(c0)

\2. c0 += random constant

\3. c2 *= random constant

\4. tmp = c1

\5. ea.shift(tmp, random amountin [-nSlots/2, nSlots/2])

\6. c2 += tmp

\7. ea.rotate(c2, random amountin [1-nSlots, nSlots-1])

\8. c1.negate()

\9. c3.multiplyBy(c2)

\10. c0 -= c3

**************/

int main(int argc, char**argv)

{

​    //long pp=101;//plaintext base

​    //long rr=1;//R-cipher times

​    //long LL=4;//levels in the modulus chain

​    //long cc=2;//number of columns in thekey-switching matrices

​    //long kk=80;//security param

​    //long ss=0;//minimum number of slots  [ default=0 ]

​    //long dd=0;// degree of the fieldextension

​    //long ww=64;//

​    //long mm = FindM(kk, LL, cc, pp, dd, ss,0);//

 

  long mm,pp,rr,LL,cc,ww,dd,kk,ss;

   ww=64;

   rr=2;

   pp=3;

   cc=2;

   mm=0;

   kk=80;

   ss=0;

   dd=1;

   LL=1+NextPowerOfTwo(dd);

   //if (mm<2)

   //mm = FindM(/*secprm=*/80, LL, /*c=*/3, pp,1, 0, mm, true);

  mm = FindM(kk, LL, cc, pp, dd, ss, 0);

 

  FHEcontext context(mm,pp,rr);//param setting

  buildModChain(context,pp,rr); //"newadding"

  cout<<context<<endl;

  //test receive

  fstream keyFile("/home/lowell/test.txt",fstream::out|fstream::trunc);//add test.txt fstream::out|fstream::trunc

  assert(keyFile.is_open());

  fstreamkeyFile1("/home/lowell/test1.txt", fstream::out|fstream::trunc);

  assert(keyFile1.is_open()); //save ciphertext

  //

  //buildModChain(context, LL, cc);

  ZZX G = context.alMod.getFactorsOverZZ()[0];// a value G which is the monic polynomial used in constructing the ring BGVworks

 

  //Key Genneration

  FHESecKey secretKey(context);

  const FHEPubKey& publicKey = secretKey;

  secretKey.GenSecKey(ww); // AHamming-weight-w secret key

  addSome1DMatrices(secretKey);

 

  //Encryption

  EncryptedArray ea(context, G);

 

  //long nslots=ea.size();

  //Use a PlaintextArray to represent aplaintext

 

  NewPlaintextArray p1(ea);//++

  random(ea,p1); //random param

  p1.print(keyFile);

 

  //cout<<p0<<endl;

  Ctxt c(publicKey);//++

  //ZZX ploy=to_ZZX(6);

  //ea.encode(ploy,p1);//++

  ea.encrypt(c, publicKey, p1);//++//plaintxtvalues p1, ciphertxt values c

  //c.print(keyFile);

  keyFile1<<c;

  //cout<<c;

  //cout<<c.getContext();// same!



  //Decryption

  NewPlaintextArray p_decrypted(ea);//++

  ea.decrypt(c, secretKey,p_decrypted);//++decrypted ciphertxt c ,nonally p_decrypted == p1

  //cout<<c;

  //cout<<c.getPtxtSpace();c isciphertext!

  //cout<<c.getContext();

 if (equals(ea,p1,p_decrypted))

​     cout<<"equal!"<<endl;

 else

​     cerr<<"oops0"<<endl;

}

###4、The following figure shows the code displayed in the experiment result.

Results of the% aerial data test

%Frist one the relationship between the dimensionalitycompression ratio and the

%reconstruction error

%Line chart

X=0:1:10;

Y=[0.2653,0.1661,0.1040,0.0739,0.0464,0.0290,0.0206,0.0141,0.0089,0.0033,0];%1，5，15，25，40，55，65，75，85，100，110

Z=[0.0057,0.0283,0.0849,0.1415,0.2263,0.3112,0.3678,0.4244,0.4809,0.5658,0.6224];

plot(X,Y,'r-v',X,Z,'b-s');

xlabel('Experiment No.');

ylabel('Ratio');

legend('The Reconstruction Error Ratio','Dimensionality Compression Ratio');

 

%comparsion between the traditional SVDdecomposition and S-tSVD method

X1=[6,600,6000,60000,600000]; %,6000000

%X = mapminmax(X1, 0, 1);

X = log10(X1);

%display(X);

Y=[0.016323,0.026745,0.027036,0.701844,5.219157];

% Y = mapminmax(Y1, 0, 1);

% display(Y);

Z=[0.010292,0.010729,0.015893,0.053003,0.353426];%,1.760767

% Z = mapminmax(Z1, 0, 1);

% display(Z);

X2=[600000,6000000];

X_2=log10(X2);

Y1=[5.219157,10];

Z1=[0.353426,1.760767];

plot(X,Y,'r-v',X,Z,'b-s');

hold on

plot(X_2,Y1,'r--',X_2,Z1,'b');

text(log10(600000),5.219157,'Exceed maximum memory');

xlabel('Logarithmic Tensor Size');

ylabel('Decomposition Time (Seconds)');

legend('The traditional SVD method','The S-tSVD method');

 

%compression ratio and RSE comparisonfor videos

 %Comparison line chart

X=0.1:0.1:0.6;

%Y=[-15.0127,-18.5740,-21.6505,-24.5485,-26.8662,-29.5708];

%Z=[-19.9176,-26.2755,-32.4904,-39.1371,-47.0240,-64.3368];

% W=[-14.5054,-18.2516,-20.9925,-23.4936,-25.3934,-27.5774];

%Y=[-16.4343,-19.2691,-21.3972,-23.3380,-25.0623,-26.9619];

%Z=[-20.6316,-25.3376,-30.1980,-35.6440,-42.3254,-56.8051];

%W=[-16.5202,-19.0675,-21.0457,-22.8575,-24.2380,-25.8724];

% Y=[-15.9127,-19.5394,-22.4385,-25.2962,-27.7576,-30.6231];

%Z=[-21.3147,-27.7269,-34.3608,-41.4481,-50.0613,-67.6347];

%W=[-16.0472,-19.8990,-23.1944,-26.0432,-28.2895,-30.8966];

Y=[-12.4470,-14.9547,-17.0013,-19.0158,-20.7456,-22.6958];

Z=[-16.7889,-21.6937,-26.5975,-32.2111,-39.0548,-54.0969];

W=[-12.4928,-14.8974,-16.8061,-18.6372,-20.0807,-21.7825];

plot(X,Y,'r-*',X,Z,'b-v',X,W,'g->');

xlabel('Compression ratio');

ylabel('RSE in db');

set(gca, 'XTick', [0.1 0.2 0.3 0.4 0.5 0.6]);

%title('Robots Battle');

%title('Aerial Cityscape');

%title('Toy Story');

title('Basketball Match');

%plot(X,Z,'-s','linewidth',2,'markersize',4);

legend('SVD compression','S-tSVD compression','tSVD-tubal compression');

 

%Illustration of four matrices of theeigenvalue tensor

obj = VideoReader('Aerial_data.mp4');

numFrames = obj.NumberOfFrames;% Total number of frames 7495

n=3; %This can be set by yourself.

D=ones(110,128,n);%Create an empty 3D matrix

i=1;

array=[]; %This array is defined to record the value of the captured frame.

for k = 1: 49: numFrames% If the number of captured frames is the same as that of n, 500: 2500 1:49

​    array(i)=k; 

​    frame = read(obj,k);%Indicates the number of the frame to be read.

Every time you read a picture, convert him into a gray picture

​    A=rgb2gray(frame); 

​    J=imresize(A,[110,128]);

​    A1=double(J);

​    D(:,:,i)=A1;

​    i=i+1;

end

x_hat=fft(D,[],3);

for i=1:size(D,n)

   [U1,S1,V1]=svd(x_hat(:,:,i));

​    if i==1

​       D1=ones(size(U1,1),size(U1,2),size(D,3));

​       D2=ones(size(S1,1),size(S1,2),size(D,3));

​       D3=ones(size(V1,1),size(V1,2),size(D,3));

​    end

   D1(:,:,i) = U1;

   D2(:,:,i) = S1;

   D3(:,:,i) = V1;

end

U=ifft(D1,[],3);

S=ifft(D2,[],3);

V=ifft(D3,[],3);   

S_matrix1=S(:,1,:);

S_matrix2=S(:,2,:);

S_matrix3=S(:,3,:);

S_matrix4=S(:,4,:);

S_vactor1=S_matrix1(:); %42240*1

S_vactor2=S_matrix2(:);

S_vactor3=S_matrix3(:);

S_vactor4=S_matrix4(:);

I_index_number1=find(S_vactor1~=0); %..*1

I_index_number2=find(S_vactor2~=0);

I_index_number3=find(S_vactor3~=0);

I_index_number4=find(S_vactor4~=0);

x=1:153;

plot(x,S_vactor1(I_index_number1,1),'.-r','MarkerSize',10);

hold on

plot(x,S_vactor2(I_index_number2,1),'.-b','MarkerSize',10);

hold on

plot(x,S_vactor3(I_index_number3,1),'.-k','MarkerSize',10);

hold on

plot(x,S_vactor4(I_index_number4,1),'.-g','MarkerSize',10);

xlabel('Experiment');

ylabel('Singular Value');

%set(gca,'xTick',1:160);

legend('Slice one','Slice two','Slice three','Slice four');

 max_S21=max(max(S(:,21,:)));

 display(max_S21);

 

 min_S12=min(min(S(:,:,12)));

 display(min_S12);



  a1=1:1:50;

   b1=1:1:50;

   [x1,y1]=meshgrid(a1,b1);

   M_S001=ones(50,50); %110*128

​    %display(S(1,1,1));

​    for j=1:50

​       for k=1:50

​       M_S001(k,j)=S(k,j,21);

​        end

​    end

   p001=mesh(x1,y1,M_S001);

   c=mean(M_S001(:)); %Average

   display(c)

   title('S(:,:,21)')