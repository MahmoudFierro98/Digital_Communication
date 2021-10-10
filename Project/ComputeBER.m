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
N_bits = 4000;
ERROR = 0;
%%% WRITE YOUR CODE HERE
for (k =1 : N_bits )
    if ( bit_seq (k)~= rec_bit_seq (k))
        ERROR = ERROR + 1 ;
    end
end 

BER = ( ERROR / N_bits ) ;
