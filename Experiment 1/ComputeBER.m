function BER = ComputeBER(bit_seq,rec_bit_seq)
%
% Inputs:
%   bit_seq:     The input bit sequence
%   rec_bit_seq: The output bit sequence
% Outputs:
%   BER:         Computed BER
%
% This function takes the input and output bit sequences and computes the
% BER

%%% WRITE YOUR CODE HERE
error_bits_number = 0;
bits_number       = length(bit_seq);
for i = 1:1:bits_number
   if  bit_seq(1,i)~= rec_bit_seq(1,i)
       error_bits_number = error_bits_number + 1;
   end
end
BER = error_bits_number / bits_number; 
%%%
