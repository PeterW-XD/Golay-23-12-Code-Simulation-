%% 子函数-格雷码解码
function C_result = decode(R,G,H,SE)
%{
输入：
    接收序列R
    生成矩阵G
    校验矩阵H
    全部伴随式与错误图样SE 
输出：
    译码结果C_result
%}
[k,n] = size(G);
S = mod(R*(H.'),2);         % 生成伴随式S
[S_row,~] = size(S);        %row 组数 column n-k

SE = mat2cell(SE,ones(1,2048),[11,23]);  %分成元胞，便于查找，相当于“查找表”

 % 找出计算出的伴随式所对应的错误图样，并进行纠正
C_result = zeros(S_row,n);          %row组数
[SE_row,~] = size(SE);              %row伴随式数2048 
 parfor m=1:S_row                   %遍历所有组数
     for n=1:SE_row                 %遍历查找表        
         if all(S(m,:) == cell2mat(SE(n,1)))    %找到对应伴随阵
             C_result(m,:) = R(m,:)+cell2mat(SE(n,2));      %输出对应错误图样+接收信号：C=R+e
             C_result(m,:) = mod(C_result(m,:),2);          %模2
         else
             C_result(m,:) = mod(C_result(m,:),2);
         end        
     end  
 end

C_result = C_result(:,1:k);                     %取前k位即信息位，即完成解码
end
