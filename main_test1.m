%%%%%%%%%%%%%%%%%%%%%%%%%%%% Golay码的编译码仿真 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% main_test1.m %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2022年05月26日 %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%程序说明
%完成了golay码在点对点通信系统中的性能分析
%设置参数如下：
%信道编码方式：(23,12)golay码  %%%%%  调制方式：BPSK
%噪声：加性高斯白噪声

%%仿真环境 matlab R2020a
%{
仿真：假如是某CDMA体制通信系统，每个用户的数据传送速率为 1kb/s，
信道编码采用编码速率为 R = 12/23 的格雷编码，采用2PSK调制方式，
基带脉冲成型采用因子为α等于 0.35 的平方根升余弦滚降函数，
上采样倍数 insValue 为 8，直序扩频倍数 d = 16 ，系统中的用户个数U = 3。
此时SNR = EbN0+10log10(12/23*1*3/16/8)
%}
%参考资料 https://blog.csdn.net/weixin_46258766/article/details/118532729
%张少侃,吕聪敏,甘浩.数字通信系统中Eb/N0与SNR转换方法的研究[J].现代计算机,2019(12):33-36.

%代码使用parfor并行计算，加快计算速度，需要先安装工具箱

%% 预计算
clear;
clc;
%  Golay(23,12,7)的P矩阵
P=[1 1 0 0 0 1 1 1 0 1 0;
   0 1 1 0 0 0 1 1 1 0 1;
   1 1 1 1 0 1 1 0 1 0 0;
   0 1 1 1 1 0 1 1 0 1 0;
   0 0 1 1 1 1 0 1 1 0 1;
   1 1 0 1 1 0 0 1 1 0 0;
   0 1 1 0 1 1 0 0 1 1 0;
   0 0 1 1 0 1 1 0 0 1 1;
   1 1 0 1 1 1 0 0 0 1 1;
   1 0 1 0 1 0 0 1 0 1 1;
   1 0 0 1 0 0 1 1 1 1 1;
   1 0 0 0 1 1 1 0 1 0 1];
% Golay(23,12,7)的生成矩阵
G = [eye(12) P];
% Golay(23,12,7)的校验矩阵
H = [transpose(P) eye(11)];

%生成所有可能的码组，用于软判决
x1=[0 1];
x2=[0 1];
x3=[0 1];
x4=[0 1];
x5=[0 1];
x6=[0 1];
x7=[0 1];
x8=[0 1];
x9=[0 1];
x10=[0 1];
x11=[0 1];
x12=[0 1];
[x12,x11,x10,x9,x8,x7,x6,x5,x4,x3,x2,x1] = ndgrid(x12,x11,x10,x9,x8,x7,x6,x5,x4,x3,x2,x1);
combin=[x1(:) x2(:) x3(:) x4(:) x5(:) x6(:) x7(:) x8(:) x9(:) x10(:) x11(:) x12(:)];            %列出所有可能12位码
combin_stream = reshape(combin.',[1,length(combin(:))]);                                        % 转换成码流

golay_combin = golay(combin_stream,G);                                                          %计算所有可能的23位golay

%H 23*11
%错一位时的伴随式和错误图样
e1 = [zeros(1,23);eye(23)]; %24*23
S1 = mod(e1*(H.'),2);  %24*11

%错两位时的伴随式和错误图样
cursor = 1;
e2 = zeros(253,23); %253,23
for j=1:23
    for k=j+1:23
        e2(cursor,j)=1;
        e2(cursor,k)=1;
        cursor=cursor+1;
    end
end
S2 = mod(e2*(H.'),2);  %253*11
clear j k;

%错三位时的伴随式和错误图样
cursor = 1;
e3 = zeros(1771,23);
for x=1:23
    for y=x+1:23
        for z=y+1:23   
        e3(cursor,x)=1;
        e3(cursor,y)=1;
        e3(cursor,z)=1;
        cursor=cursor+1;
        end
    end
end
S3 = mod(e3*(H.'),2);  
clear x y z cursor;

%最终的伴随式与错误图样  
SE = [S1 e1;S2 e2;S3 e3];


%% 主函数

N = input('每个信噪比条件下要发送多少个样本点：');

NbN0_min = 0; NbN0_max = 11;% NbN0范围是0-11dB
NbN0 = NbN0_min:0.5:NbN0_max;
snrdB = (NbN0_min:0.5:NbN0_max)+10*log10(((12/23)*3)/16*8); %转换为SNR
sym_initial = round(rand(1,N));                             %生成原始序列
% 获得格雷码编码C_golay
[k,n] = size(G);
C_golay = golay(sym_initial,G);                             %用golay函数生成格雷码
[groups,~] = size(C_golay);                                 %检查线性分组码分组数
C_golay_stream = reshape(C_golay.',[1,length(C_golay(:))]); %转换成码流


tic
% 两种序列都通过BPSK-AWGN信道
[~,BER_dir] = BPSK(sym_initial,snrdB);          % 未编码码流
[RX_golay,~] = BPSK(C_golay_stream,snrdB);      % 格雷编码码流（硬判决）
RX_golay_s = BPSK_soft(C_golay_stream,snrdB);   % 格雷编码码流（软判决）输出
BER_golay_s = BER_golay_soft(RX_golay_s,combin,golay_combin,C_golay,snrdB,n,groups,N);% 格雷编码码流（软判决）误码率
toc

tic
% 格雷码解码并计算SNR
errors_golay = zeros(1,length(snrdB));  % 预分配错误内存
for i=1:length(snrdB)
    RX_golay1 = reshape(RX_golay(i,:),[n,groups]).';                    %重整
    
    C_result1 = decode(RX_golay1,G,H,SE);                               %解码---------------high delay
    
    errors_golay(i) = sum(sum(mod(C_result1+C_golay(:,1:k),2)));        %计算误码数
end
BER_golay = errors_golay/N;                                             %得到BER
toc

% 画出三种方案的BER-SNR曲线
figure
semilogy(NbN0,BER_dir,'*-',NbN0,BER_golay,'o-',NbN0,BER_golay_s,'+-');grid on;
axis([NbN0_min NbN0_max 1e-6 1]);
xlabel('NbN0/ dB');ylabel('误码率 BER');
title(['发送的信息序列长度为 ',num2str(N)]);
legend('未编码系统','格雷编码系统（硬判决）','格雷编码系统（软判决）');
% 直接显示数值
disp(['未编码系统',num2str(BER_dir)]);
disp(['格雷编码系统',num2str(BER_golay)]);
disp(['格雷编码系统（软判决）',num2str(BER_golay_s)]);

%% 最小欧氏距离准则解码加误码率计算
function BER = BER_golay_soft(RX,combin1,combin2,C_golay,SNR,n,groups,N)
%{
输入：
    接收信号RX
    标准原始分组矩阵 查找表1 combin1
    标准编码后矩阵 查找表2 combin2
    编码后原始矩阵C_golay
    信噪比SNR
    格雷码长n
    组数groups
    总bit数N
输出：
    误码率BER
%}
errors = zeros(1,length(SNR));
result = zeros(groups,12);
combin2 = 2*combin2 - 1;
for i=1:length(SNR)
    RX_golay1 = reshape(RX(i,:),[n,groups]).';              %重整
    parfor j=1:groups
        R_soft = RX_golay1(j,:);
        R_soft1 = repmat(R_soft,4096,1);                    %生成每行都相同的，减少循环体
        d = sqrt(sum((R_soft1 - combin2).^2,2));            %每行加起来
        [~,num] = min(d);                                   %找出最小的距离位置
        result(j,:) = combin1(num,:);                       %读出最小距离对应的码
    end 
    errors(i) = sum(sum(mod(result+C_golay(:,1:12),2)));    %统计误码数
end
BER = errors/N;                                             %得到BER 
end


%% 子函数-完成BPSK-AWGN信道的仿真（软输出）
function RX = BPSK_soft(TX,SNR)
%{
输入
    带发送序列TX
    信噪比SNR
输出
    发送信号（未判决）RX
%}
N = length(TX);             % 获得原始序列长度
len_snr = length(SNR);      % 获得SNR范围长度

RX = zeros(len_snr,N);      % 预分配接收判别序列内存

for j=1:len_snr             % 遍历所有SNR
    x = SNR(j);
   parfor i=1:N                        % 遍历每一个发送符号
      x_d = 2*TX(i) - 1;            % 得到+1和-1的发送序列
      y_d = awgn(x_d,x ,'measured');	% 加性噪声
      RX(j,i) = y_d;                        %soft
   end
end
end

%% 其他函数在文件夹中