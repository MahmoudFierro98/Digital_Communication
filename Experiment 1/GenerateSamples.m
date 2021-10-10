function sample_seq = GenerateSamples(bit_seq,fs)
%
% Inputs:
%   bit_seq:    Input bit sequence
%   fs:         Number of samples per bit
% Outputs:
%   sample_seq: The resultant sequence of samples
%
% This function takes a sequence of bits and generates a sequence of
% samples as per the input number of samples per bit

sample_seq = zeros(size(bit_seq*fs));

%%% WRITE YOUR CODE FOR PART 2 HERE
sample_seq = zeros(1,length(bit_seq)*fs);
n = 1;
for i = 1:1:(length(bit_seq)*fs)
    sample_seq(1,i) = bit_seq(1,n);
    if mod(i,fs) == 0
        n = n + 1;
    end
end
%%%