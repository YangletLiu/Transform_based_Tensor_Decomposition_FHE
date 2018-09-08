%代码的整合阶段
%我们的实验数据是video。有关Aerial data（航拍数据）、感知数据sensory data （Intelligent building和Intelligent Transportation case）
%涉及全同态加解密实验 
%这个首先要在linux配置中设置好HElib库，然后测试Test_General的代码

%涉及的基本代码
obj = VideoReader('Aerial_data.mp4');
numFrames = obj.NumberOfFrames;% 帧的总数 7495
n=3; %这个可以自己设定
D=ones(110,128,n);%创建一个空的三维矩阵
i=1;
array=[]; %定义这个数组就是为了记录截取的帧的值
for k = 500 : 2500 : numFrames% 保持截出来的帧的数量和n保持一致就行 500 : 2500 1:49
     array(i)=k; 
     frame = read(obj,k);%读取第几帧
     %每读一张把他转化成灰度图片
     A=rgb2gray(frame); 
     J=imresize(A,[110,128]);
     A1=double(J);
     D(:,:,i)=A1;
     i=i+1;
     %figure(k);
     %imshow(A);
     %-----------
     %figure(k);
     %imshow(frame);%显示帧
     %imwrite(J,strcat('/Users/apple/Desktop/1/',num2str(k),'.jpg'),'jpg');% 保存帧
end
%display(array);
B1=[D(:,:,2);D(:,:,3);D(:,:,1)];
%imshow(uint8(B1));
%从文件中读取你想选择的照片---读取照片（一）------- 单个的一张照片
% [filename, pathname] = uigetfile('*.jpg', '读取图片文件'); %选择图片文件
% if isequal(filename,0)   %判断是否选择
%    msgbox('没有选择任何图片');
% else
%    pathfile=fullfile(pathname, filename);  %获得图片路径
%    M=imread(pathfile);     %将图片读入矩阵
%    %imshow(M);    %绘制图片
% end
%循环读入一个文件中所有的照片-----读取照片（二）------- M存放了所有的整个视频的所有照片
% fileform = '/Users/apple/Desktop/1/*.jpg';
% filepathsrc = '/Users/apple/Desktop/1/';
% file = dir(fileform); %获得所有指定文件夹下的所有子文件夹和文件，并存放在一种为文件结构体数组中
% M=[];
% for i = 1:length(file)
%         m = imread([filepathsrc, file(i).name]);
%         M = [M;m];
% end
%----------------------------------------------------------------------------
%T-SVD分解
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
%计算compression degree(Dimensionality Reduction Ratio)
%---估算逼近程度（降维实验）
    %选取前面__个非零奇异值  原始110*128
    for n1 = 1:5  %（1）
        for n2 = 1:5 %（2）
            for n3=1:3
            S_1(n1,n2,n3) = S(n1,n2,n3);  
            end
        end  
    end 
    %左奇异张量
    for n1 = 1:110  
        for n2 = 1:5  %（3）
            for n3 = 1:3
            U_1(n1,n2,n3) = U(n1,n2,n3);  
            end
        end  
    end  
    %右奇异张量
    for n1 = 1:128
        for n2 = 1:5 %（4）
            for n3= 1:3
            V_1(n1,n2,n3) = V(n1,n2,n3);  
            end
        end  
    end  
   %计算张量的product??New Tensor Multiplication
   %display(size(S_1));
    c = t_product(U_1,S_1); %c计算出来是个tensor
    % V_1求转置
    h2 = ones(5,128,3); %（5）-1位置
    for j=1:size(D,n)
        h2(:,:,j)=V_1(:,:,j)'; %110*126*3
    end
    for k=2:size(D,n)
        t1=h2(:,:,k);
        h2(:,:,k)=h2(:,:,2+2-k);
        h2(:,:,2+2-k)=t1;
    end
   c = t_product(c,h2);
   result=[c(:,:,1);c(:,:,2);c(:,:,3)];
   %imshow(uint8(result));
   %---------
    %求相对误差
    T_error=result-B1;
    display(T_error);
    %计算重构误差率
    n = norm(T_error, 'fro' );
    %display(n);
    %原始的张量的Frobenius范数
    n_B=norm(B1,'fro');
    %Tensor Approximation Ratio(reconstruction error degree)
    e=n/n_B; 
    display(e);
    RES=20*log10(e);
    display(RES);
    p = size(D,1)*size(D,2)*size(D,3)/(n2*(size(D,1)+size(D,2)+1));
    display(1/p);
%-------------------------------------
%涉及和traditional svd、tSVD-slice分解的比较
    %traditional svd method
    addpath /Users/apple/Downloads/tensor_toolbox 
    obj = VideoReader('Aerial_data.mp4');
    numFrames = obj.NumberOfFrames;% 7495
    D=ones(110,128,3); %创建一个空的三维矩阵
    i=1;
    for k = 500 : 2500 : numFrames
         frame = read(obj,k);%读取第几帧
         %每读一张把他转化成灰度图片
         A=rgb2gray(frame); 
         J=imresize(A,[110,128]);
         A1=double(J);
         D(:,:,i)=A1;
         i=i+1;
    end
    A=tensor(D);
    tic
    A1=tenmat(A,1);  %沿着mode-1展开成矩阵 n1*(n2*n3)
    %A1=A1';
    [U,S,V]=svd(A1.data); 
    %toc
    %选取前面__个非零奇异值，乘起来
    %tic
    U1=U(:,1:53); %(1)修改
    S1=S(1:53,1:53); %(2)修改
    V1=V(:,1:53); %(3)修改
    AA=U1*S1'*V1';
    toc
    A1=double(A1);  
    SVD_Error = AA-A1;
    n_svd = norm(SVD_Error, 'fro' );
    n_B_svd = norm(A1,'fro');
    e_svd = n_svd/n_B_svd; 
    display(e_svd);
    RES_svd=20*log10(e_svd);
    display(RES_svd);
    %--------------------
    %tSVD-tubal method
    obj = VideoReader('Aerial_data.mp4');
    numFrames = obj.NumberOfFrames;
    n=3;
    D=ones(110,128,n);%创建一个空的三维矩阵
    i=1;
    for k = 500 : 2500 : numFrames
         frame = read(obj,k);%读取第几帧
         A=rgb2gray(frame); 
         J=imresize(A,[110,128]);
         A1=double(J);
         D(:,:,i)=A1;
         i=i+1;
    end
    tic %无fft频域变化  t-product算法
     for i=1:size(D,n)
        [U1,S1,V1]=svd(D(:,:,i));
        if i==1
            D1=ones(size(U1,1),size(U1,2),size(D,n));
            D2=ones(size(S1,1),size(S1,2),size(D,n));
            D3=ones(size(V1,1),size(V1,2),size(D,n));
        end
        D1(:,:,i) = U1;
        D2(:,:,i) = S1;
        D3(:,:,i) = V1;
     end
     %实行降维处理
     for n1 = 1:35 %（1）
        for n2 = 1:35 %（2）
            for n3=1:3
            S_1(n1,n2,n3) = D2(n1,n2,n3);  
            end
        end  
    end 
    %左奇异张量
    for n1 = 1:110  
        for n2 = 1:35 %（3）
            for n3 = 1:3
            U_1(n1,n2,n3) = D1(n1,n2,n3); 
            end
        end  
    end  
    %右奇异张量
    for n1 = 1:128
        for n2 = 1:35 %（4）
            for n3=1:3
            V_1(n1,n2,n3) = D3(n1,n2,n3);
            end
        end  
    end  
%     %矩阵乘法
    for m = 1:1:3
        mid(:,:,m) = U_1(:,:,m)*S_1(:,:,m);
    end
    for n=1:1:3
        result1(:,:,n)= mid(:,:,n)*V_1(:,:,n)';
    end
    toc
    %new tensor product
%     c = t_product(U_1,S_1); %c计算出来是个tensor
%     % V_1求转置
%     h2 = ones(35,128,3); %（5）-1位置
%     for j=1:size(D,n)
%         h2(:,:,j)=V_1(:,:,j)'; %110*126*3
%     end
%     for k=2:size(D,n)
%         t1=h2(:,:,k);
%         h2(:,:,k)=h2(:,:,2+2-k);
%         h2(:,:,2+2-k)=t1;
%     end
%    c = t_product(c,h2);
%    result=[c(:,:,1);c(:,:,2);c(:,:,3)];
    %-------结果的呈现----------
    B1=[D(:,:,2);D(:,:,3);D(:,:,1)];
    %imshow(uint8(B1));
    %imshow(uint8(D));
    result=[result1(:,:,2);result1(:,:,3);result1(:,:,1)];
    %imshow(uint8(result));
    %求相对误差
    T_error=result-B1;
    %display(T_error);
    %计算重构误差率
    n = norm(T_error, 'fro' );
    n_B=norm(B1,'fro');
    e=n/n_B; 
    display(e);
    RES=20*log10(e);
    display(RES);
%----------------------------------------------
%涉及计算运行时间running time compare traditional svd
    %traditional svd time
    addpath /Users/apple/Downloads/tensor_toolbox 
    obj = VideoReader('Aerial_data.mp4');
    numFrames = obj.NumberOfFrames; 
    D=ones(100,200,30);%创建一个空的三维矩阵 %第一处要修改的地方
    i=1;
    for k = 500 : 2500 : numFrames%
         frame = read(obj,k);%读取第几帧
         %每读一张把他转化成灰度图片
         A=rgb2gray(frame); 
         J=imresize(A,[100,200]); %第二处要修改的地方
         A1=double(J);
         D(:,:,i)=A1;
         i=i+1;
    end
    A=tensor(D);
    tic
    A1=tenmat(A,1);  
    [U,S,V]=svd(A1.data); 
    toc  %[1,2,3]:0.016323  [10,20,3]:0.026745 [10,20,30]:0.027036  [10,200,30]:0.701844 [100,200,30]:5.219157  [100,200,300]Error using svd
    %Requested 60000x60000 (26.8GB) array exceeds maximum array size preference. Creation of arrays greater than
    %this limit may take a long time and cause MATLAB to become unresponsive. See array size limit or preference
    %panel for more information.
    %S-tSVD
    obj = VideoReader('Aerial_data.mp4');
    numFrames = obj.NumberOfFrames;% 9243 
    D=ones(10,200,30);%创建一个空的三维矩阵
    i=1;
    for k = 1 : 309: numFrames% 3081 309 31
         frame = read(obj,k);%读取第几帧
         %每读一张把他转化成灰度图片
         A=rgb2gray(frame); 
         J=imresize(A,[10,200]);
         A1=double(J);
         D(:,:,i)=A1;
         i=i+1;
    end
    x_hat=fft(D,[],3); %任意循环矩阵可以被傅里叶变换矩阵对角化
    tic
    for i=1:size(D,3)
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
    toc %[1,2,3]:0.010292 [10,20,3]:0.010729 [10,20,30]:0.015893 [100,200,30]:0.353426 [100,200,300]:1.760767 [1000,2000,3000]:22.595927
    U=ifft(D1,[],3);
    S=ifft(D2,[],3);
    V=ifft(D3,[],3);
%------------------------------------------
%涉及求得张量分解之后的特征矩矩阵、特征值
    max_S2=max(max(S(:,:,2)));
    display(max_S2);
    min_S22=min(min(S(:,:,2)));
    display(min_S22);
    %S（：，；，k）中每一个切片的信息
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
    %S（：，k，：）的切片信息
     T_S2=S(:,2,:);
     display(T_S2); %110*1*155
     M_S2=ones(110,155);
    %  for i=1:155
    %       M_S2(:,i)=T_S2(:,:,i); %110行1列
    %  end
     for i=1:110
         for j=1:155
             M_S2(i,j)=T_S2(i,1,j);
         end
     end
     c2=mean(M_S2(:)); %全部平均
     display(c2); 
     p2=mesh(x,y,M_S2);
     title('S(:,2,:)');
     
     %S（：，j，：）的切片信息
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
     x=1:155;
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

