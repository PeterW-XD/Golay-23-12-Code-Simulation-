%% 子函数-格雷码编码
function C = golay(M,G)
%{
输入：
    发送序列M
    生成矩阵G
输出：
    编码结果C
%}
[k,~] = size(G);             %行数

% 输入序列补位
N = size(M,2);               % 获得输入序列元素个数
r = mod((k-rem(N,k)),k);     % 获得需要对输入序列进行补位的个数
M_add0 = [M,zeros(1,r)];     % 补位

% 将输入信息序列进行分组
groups = ceil(length(M_add0)/k);        %获得分组个数
M_dis = reshape(M_add0,[k,groups]).';   %重整

% 生成编码结果C
C = mod(M_dis*G,2);                     %通过生成矩阵生成格雷码，需要模2
end
