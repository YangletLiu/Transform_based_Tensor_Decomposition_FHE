%航拍数据实验结果显示
%第一个The relationship between the dimensionality compression ratio and the
%reconstruction error
%折线图
X=0:1:10;
Y=[0.2653,0.1661,0.1040,0.0739,0.0464,0.0290,0.0206,0.0141,0.0089,0.0033,0]; %1，5，15，25，40，55，65，75，85，100，110
Z=[0.0057,0.0283,0.0849,0.1415,0.2263,0.3112,0.3678,0.4244,0.4809,0.5658,0.6224];
plot(X,Y,'r-v',X,Z,'b-s');
xlabel('Experiment No.');
ylabel('Ratio');
legend('The Reconstruction Error Ratio','Dimensionality Compression Ratio');

%comparsion between the traditional SVD decomposition and S-tSVD method
X1=[6,600,6000,60000,600000]; %,6000000
%X = mapminmax(X1, 0, 1);
X = log10(X1);
%display(X);
Y=[0.016323,0.026745,0.027036,0.701844,5.219157]; 
% Y = mapminmax(Y1, 0, 1);
% display(Y);
Z=[0.010292,0.010729,0.015893,0.053003,0.353426]; %,1.760767
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

%compression ratio and RSE comparison for videos
 %对比折线图
X=0.1:0.1:0.6;
% Y=[-15.0127,-18.5740,-21.6505,-24.5485,-26.8662,-29.5708];
% Z=[-19.9176,-26.2755,-32.4904,-39.1371,-47.0240,-64.3368];
% W=[-14.5054,-18.2516,-20.9925,-23.4936,-25.3934,-27.5774];
% Y=[-16.4343,-19.2691,-21.3972,-23.3380,-25.0623,-26.9619];
% Z=[-20.6316,-25.3376,-30.1980,-35.6440,-42.3254,-56.8051];
% W=[-16.5202,-19.0675,-21.0457,-22.8575,-24.2380,-25.8724];
% Y=[-15.9127,-19.5394,-22.4385,-25.2962,-27.7576,-30.6231];
% Z=[-21.3147,-27.7269,-34.3608,-41.4481,-50.0613,-67.6347];
% W=[-16.0472,-19.8990,-23.1944,-26.0432,-28.2895,-30.8966];
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

%Illustration of four matrices of the eigenvalue tensor
obj = VideoReader('Aerial_data.mp4');
numFrames = obj.NumberOfFrames;% 帧的总数 7495
n=3; %这个可以自己设定
D=ones(110,128,n);%创建一个空的三维矩阵
i=1;
array=[]; %定义这个数组就是为了记录截取的帧的值
for k = 1 : 49 : numFrames% 保持截出来的帧的数量和n保持一致就行 500 : 2500 1:49
     array(i)=k; 
     frame = read(obj,k);%读取第几帧
     %每读一张把他转化成灰度图片
     A=rgb2gray(frame); 
     J=imresize(A,[110,128]);
     A1=double(J);
     D(:,:,i)=A1;
     i=i+1;
end
x_hat=fft(D,[],3);
for i=1:size(D,n)
    [U1,S1,V1]=svd(x_hat(:,:,i));
    if i==1
        D1=ones(size(U1,1),size(U1,2),size(D,3));
        D2=ones(size(S1,1),size(S1,2),size(D,3));
        D3=ones(size(V1,1),size(V1,2),size(D,3));
    end
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
%
   a1=1:1:50;
    b1=1:1:50;
    [x1,y1]=meshgrid(a1,b1);
    M_S001=ones(50,50); %110*128
    %display(S(1,1,1));
    for j=1:50
        for k=1:50
        M_S001(k,j)=S(k,j,21);
        end
    end
    p001=mesh(x1,y1,M_S001);
    c=mean(M_S001(:)); %全部平均
    display(c)
    title('S(:,:,21)')





