%子函数-完成BPSK-AWGN信道的仿真
function [RX,BER] = BPSK(TX,SNR)
%{
输入：
    原始码元序列Tx
    信噪比范围SNR(dB)
输出：
    接收判别后序列RX
    信噪比对应的误码率序列BER
%}

N = length(TX);             % 获得原始序列长度
snr = 10.^(SNR/10);         % 转化成公制
len_snr = length(snr);      % 获得SNR范围长度

RX = zeros(len_snr,N);      % 预分配接收判别序列内存
errors = zeros(1,len_snr);  % 预分配错误内存
for j=1:len_snr             % 遍历所有SNR
   error_count = 0;         %error计数置零
   x = SNR(j);          
   parfor i=1:N                         % 遍历每一个发送符号
      x_d = 2*TX(i) - 1;                % 得到+1和-1的发送序列
      y_d = awgn(x_d,x,'measured');     % AWGN信道
      if y_d > 0
         RX(j,i) = 1;
      else
         RX(j,i) = 0;
      end
      if (RX(j,i) ~= TX(i))
         error_count = error_count + 1;	% 对错误样本进行计数
      end
   end
   errors(j) = error_count;             % 得到该信噪比下的错误个数
end
BER = errors/N;                         % BER计算
end
