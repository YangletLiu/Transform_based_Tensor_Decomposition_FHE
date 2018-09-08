%t_product函数实现过程
function c=t_product(a,b)
%     a_n1=110;a_n2=110;a_n3=3;
%     b_n1=110;b_n2=110;b_n3=3;
    [a_n1 a_n2 a_n3]=size(a); %n1 n2 n3 values 110 110 3
    [b_n1 b_n2 b_n3]=size(b); %110 110 3
    c=zeros(a_n1,b_n2,a_n3); 
    A=cell(a_n3,1);
    B=cell(b_n3,1);
    
    for i=1:a_n3
        A{i}=a(:,:,i);
        B{i}=b(:,:,i);
    end
    index_up=zeros(1,a_n3);
    index_down=zeros(1,a_n3);
    for i=1:a_n3    
        index_up(i)=a_n3-i+1;
        index_down(i)=i;
    end
    %s存放了展开的tensor的左项
    s=cell(a_n3,a_n3);
    for i=1:a_n3
        for j=1:a_n3
            if i==j
                s{i,j}=A{1};
            end       
            if j>i
                s{i,j}=A{index_up(j-i)};
            end       
            if j<i
                s{i,j}=A{index_down(i-j+1)};
            end      
        end   
    end
    %申请一个存放计算结果的内存大小
    re=cell(a_n3,1);
    for i=1:a_n3
        re{i}=zeros(a_n1,b_n2);
    end
%     display(re);

    for i=1:a_n3
        for j=1:a_n3
            for k=1:1
                re{i,k}=re{i,k}+s{i,j}*B{j,k};
            end
        end    
    end
% 计算求得乘积的解
    for i=1:a_n3
        c(:,:,i)=re{i};        
    end
end