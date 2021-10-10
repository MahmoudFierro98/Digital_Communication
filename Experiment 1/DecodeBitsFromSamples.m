function rec_bit_seq = DecodeBitsFromSamples(rec_sample_seq,case_type,fs)
%
% Inputs:
%   rec_sample_seq: The input sample sequence to the channel
%   case_type:      The sampling frequency used to generate the sample sequence
%   fs:             The bit flipping probability
% Outputs:
%   rec_sample_seq: The sequence of sample sequence after passing through the channel
%
% This function takes the sample sequence after passing through the
% channel, and decodes from it the sequence of bits based on the considered
% case and the sampling frequence

if (nargin <= 2)
    fs = 1;
end

switch case_type
    
    case 'part_1'
        %%% WRITE YOUR CODE FOR PART 1 HERE
        rec_bit_seq = zeros(1,length(rec_sample_seq));
        rec_bit_seq = rec_sample_seq;
        %%%
    case 'part_2'
        %%% WRITE YOUR CODE FOR PART 2 HERE
        rec_bit_seq = zeros(1,length(rec_sample_seq)/fs);
        n = 1;
        number_of_ones  = 0;
        number_of_zeros = 0;
        for i = 1:1:length(rec_sample_seq)
            if (rec_sample_seq(1,i) == 1)
                number_of_ones = number_of_ones + 1;
            else
                number_of_zeros = number_of_zeros + 1;
            end
            if (mod(i,fs) == 0)
                if number_of_ones > number_of_zeros
                    rec_bit_seq(1,n) = 1;
                else
                    rec_bit_seq(1,n) = 0;
                end
                n = n + 1;
                number_of_ones  = 0;
                number_of_zeros = 0;
            end
        end
        %%%
    case 'part_3'
        %%% WRITE YOUR CODE FOR PART 3 HERE
        rec_bit_seq = zeros(1,length(rec_sample_seq)/fs);
        n = 1;
        for i = 1:1:length(rec_sample_seq)
            if (mod(i,fs) == 0)
                rec_bit_seq(1,n) = rec_sample_seq(1,i);
                n = n + 1;
            end
        end
        %%%
end