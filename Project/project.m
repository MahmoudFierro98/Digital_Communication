clear all;
clc;
%%
%% part-1  Inter-Symbol Interference due to band-limited channels
%%
%%%%%%%%%%%  first square input %%%%%%%%%%%%
N = 1024000-1 ;                 % Total number of samples
fs=10000000;
Ts = 1/fs;
T_Bit=2/100000;
t_axis = (0:N-1)*Ts;  
f_axis = -fs/2:fs/N:fs/2-1/N;
N_sq = round(T_Bit/Ts);
x_first = zeros(1,length(t_axis));
x_first(1:N_sq) = 1;
%x_first(1:N_sq)=ones(1,N_sq);
%x_first(N_sq+1:2*N_sq)=zeros(1,N_sq);
%t_axis=linspace(0,4e-5,N_sq*2);
figure;
plot(t_axis,x_first,'linewidth',3);
xlim([0 2*T_Bit])
ylim([0 2])
title('first square pulse before the channel in time')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%%%%%%%%%%%  first square input in freq domine %%%%%%%%%%%%
x_first_freq = fftshift(fft(x_first));
x_fir_freq_abs=abs(x_first_freq)
figure;
plot(f_axis,x_fir_freq_abs,'linewidth',3)
title('first square pulse before the channel in frequency','linewidth',10)
grid on
xlim([-100000 100000]*2)
ylim([0 200])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%%%%%%%%%%%  second square input %%%%%%%%%%%%
x_sec = ones(1,length(t_axis));
x_sec(1:N_sq) = 0;
figure;
plot(t_axis,x_sec,'linewidth',3);
xlim([0 2*T_Bit])
ylim([0 2])
title('second square pulse before the channel in time')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%%%%%%%%%%%  second square input in freq domine %%%%%%%%%%%%
x_second_freq = fftshift(fft(x_sec));
x_sec_freq_abs=abs(x_second_freq);
figure;
plot(f_axis,x_sec_freq_abs,'linewidth',3)
title('second square pulse before the channel in frequency','linewidth',10)
grid on
xlim([-100000 100000]*2)
ylim([0 200])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%channel%%%%%%%%%%%%%%%%%%%%
L=length(f_axis);
fm=100000;
n=fs/(L-1);                                                                                          
f1=-fs/2:n:-fm;                                 
len1=length(f1);                                  
yz1=zeros(1,len1);                            
f3=fm:n:fs/2;
len3=length(f3);
yz3=zeros(1,len3-1);
f2=-fm:n:fm;
len2=length(f2);
yz2=ones(1,len2);
channel=[yz1 yz2 yz3];
figure;
subplot(2,1,1);                   
plot(f_axis,channel,'linewidth',3);
title('channel','linewidth',10);
xlim([-100000 100000]*2)
ylim([0 2])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%%%%%%%%%%%%%%%%%%%%%%%channel in time %%%%%%%%%%%%%%%%%%%%
channel_time = ifft(ifftshift(channel));
figure;
plot(t_axis,channel_time,'linewidth',3)
title('channel in time','linewidth',10)
xlim([0 2*T_Bit])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


%%%%%%%%%%%%%%%%%%%%%%%first one  multipling in freq %%%%%%%%%%%%%%%%%%%%
receive_first_freq=channel.*x_first_freq;
figure;
plot(f_axis,abs((receive_first_freq)),'linewidth',3)
title(' first signal after the channel in freq')
xlim([-fm fm]*2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%%%%%%%%%%%%%%%%%%%%%%%first one  in time %%%%%%%%%%%%%%%%%%%%
receive_time = ifft(ifftshift(receive_first_freq));
figure;
plot(t_axis,receive_time,'linewidth',3)
title(' first signal after the channel in time')
xlim([0 2*T_Bit])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%%%%%%%%%%%%%%%%%%%%%%%second one  multipling in freq %%%%%%%%%%%%%%%%%%%%
receive_second_freq=channel.*x_second_freq;
figure,
plot(f_axis,abs((receive_second_freq)),'linewidth',3)
title('second signal after the channel in freq')
xlim([-fm fm]*2)
ylim([0 200])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%%%%%%%%%%%%%%%%%%%%%%% second one in time %%%%%%%%%%%%%%%%%%%%
receive_time2 = ifft(ifftshift(receive_second_freq));
figure;
plot(t_axis,receive_time2,'linewidth',3)
title(' second signal after the channel in time')
xlim([0 2*T_Bit])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%%%%%%%%%%%%%%%%%%%%%%% plotting the first and second pulses in time %%%%%%%%%%%%%%%%%%%%
figure;
plot(t_axis,x_first,'b',t_axis,x_sec,'r','linewidth',3)
title(' first and second square pulses in time')
ylim([0 2])
xlim([0 2*T_Bit])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%%%%%%%%%%%%%%%%%%%%%%% plotting the first and second pulses in time after the channel %%%%%%%%%%%%%%%%%%%%
figure;
plot(t_axis,receive_time,'b',t_axis,receive_time2,'r','linewidth',3)
title(' first and second square pulses  after the channel in time')
ylim([0 2])
xlim([0 2*T_Bit])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


%%%%%%%%%%%%%%%%%%%%%%Raised-Cosine Pulse Shaping%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%Raised-Cosine Pulse Shaping in freq %%%%%%%%%%%%%%%%%%%%%%%%
%for beta=0:.5:1
beta = 1;              % The Roll-Off Factor.
raise_cos_freq = zeros(1,N);   
ind = 1;
for i=f_axis
    if (abs(i) > (1-beta)/(2*T_Bit)) && (abs(i) <= (1+beta)/(2*T_Bit))
        raise_cos_freq(ind) = 0.5 * (1+cos((pi*T_Bit/beta)*(abs(i)-((1-beta)/(2*T_Bit)))));
        ind = ind +1;
    elseif (abs(i) <= (1-beta)/(2*T_Bit))
        raise_cos_freq(ind) = 1;
        ind = ind +1;
    else
        raise_cos_freq(ind) = 0;
        ind = ind + 1;
    end
    
end
  figure;
  plot(f_axis,raise_cos_freq,'linewidth',3);
  title('raised cos in freq')
  xlim([-100000 100000])
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%

 %%%%%%%%%%%%%%%%%%%%%%%Raised-Cosine Pulse Shaping in time %%%%%%%%%%%%%%%%%%%%%%%%
  raisedcos_time = ifft(ifftshift(raise_cos_freq));
  figure;
  plot(t_axis,raisedcos_time,'linewidth',3);
  title('raised cos in time')
  xlim([0 5*T_Bit]*2)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
  
 %%%%%%%%%%%%%%%%%%%%%% generating delta for shifiting the signal%%%%%%%%%%%%%%%%%%%%%%%%
  in1 = zeros(1,N);
  in1(1) = 1;
  in2 = zeros(1,N);
  in2(201) = 1;
  figure;
  plot(t_axis,in1,'linewidth',3);
  hold on;
  plot(t_axis,in2,'linewidth',3);
  title('delta in time')
  xlim([0 5*T_Bit])
  ylim([-1.5 1.5])  
  inf1 = fftshift(fft(in1));
  inf2 = fftshift(fft(in2));
  figure;
  plot(f_axis,inf1);
  hold on;
  plot(f_axis,inf2,'linewidth',3);
  title('delta in freq')
  xlim([-100000 100000]*2)
  ylim([-1.5 1.5])  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%shifting%%%%%%%%%%%%%%%%%%%%%%%%
  in1_r_freq = 200*raise_cos_freq .* inf1;   
  in2_r_freq = 200*raise_cos_freq .* inf2;
  in1_r1 = ifft(ifftshift(in1_r_freq));
  in2_r2 = ifft(ifftshift(in2_r_freq));
  figure;
  plot(t_axis,in1_r1,'linewidth',3);
  hold on;
  plot(t_axis,in2_r2,'linewidth',3);
  title('input to channel in time')
  xlim([0 5*T_Bit])
  figure;
  plot(f_axis,abs(in1_r_freq),'linewidth',3);
  hold on;
  plot(f_axis,abs(in2_r_freq),'linewidth',3);
  title('input to channel in frq')
  xlim([-100000 100000]*2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%% passing the channel%%%%%%%%%%%%%%%%%%%%%%
  r1_out_channel = channel .* in1_r_freq;
  r2_out_channel = channel .* in2_r_freq;
  r_1_out_time = ifft(ifftshift(r1_out_channel));
  r_2_out_time = ifft(ifftshift(r2_out_channel));
  figure;
  plot(t_axis,r_1_out_time,'linewidth',3);
  hold on;
  plot(t_axis,r_2_out_time,'linewidth',3);
  title('out of channel in time')
  xlim([0  5*T_Bit])
  ylim([-1.5 1.5])
  figure;
  plot(f_axis,abs(r1_out_channel),'linewidth',3);
  hold on;
  plot(f_axis,abs(r2_out_channel),'linewidth',3);
  title('out of channel in freq')
  xlim([-100000 100000])
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
  N_bits = 5119;
  random_Bits = randi([0 1],1,N_bits,'int8');
  deltas = zeros(1,N);
  for i=1:N_bits
      deltas(1,((i-1)*N_sq)+1) = random_Bits(1,i);
  end
  figure;
  plot(t_axis,deltas);
  title('random bits generation')
  %xlim([0 T_Bit*N_bits])
  xlim([0 T_Bit*50])
  ylim([-1 1])
%%
  deltas_f = fftshift(fft(deltas));
  %plot(f_axis,deltas_f);
  deltas_r_freq = 200*raise_cos_freq .* deltas_f;   
  deltas_r_time = ifft(ifftshift(deltas_r_freq));
  %figure;
  %plot(t_axis,deltas_r_time,'linewidth',1);
  %title('after pulse shaper in time')
  %hold on;
  %xlim([0 T_Bit*N_bits])
  %xlim([0 T_Bit*5])
  %figure;
  %plot(f_axis,deltas_r_freq,'linewidth',1);
  %title('after pulse shaper in frequency')
  %xlim([-500 500])
  % ylim([0 100000])


%%
  deltas_out_channel_freq = channel .* deltas_r_freq;
  deltas_out_channel_time = ifft(ifftshift(deltas_out_channel_freq));
%  figure;
 % plot(t_axis,deltas_out_channel_time);
  %title('out of deltas channel in time')
  %xlim([0  5*T_Bit])
  %ylim([-1.5 1.5])
 % figure;
 % plot(f_axis,deltas_out_channel_freq);
  %title('out of deltas channel in freq')
  %xlim([-100000 100000])
   %  ylim([0 100000])

%%
  %%%%%%%%%%%%%%%%%%%%%%% threshold %%%%%%%%%%%%%%%%%%%%%%
 
  
  exp_bit = zeros(1,N);
  epsilon = 0.5;
  for i=1:N_bits
      if deltas_out_channel_time(1,((i-1)*N_sq)+1) > epsilon
          exp_bit(1,i) = 1;
      else
          exp_bit(1,i) = 0;
     
      end
  end
%%
  %%%%%%%%%%%%%%%%%%%%%%% compute BER %%%%%%%%%%%%%%%%%%%%%%
  
  e = 0;
  for i=1:N_bits
      if exp_bit(1,i) ~= random_Bits(1,i)
          e = e+1;
      end
  end
  BER = e/N_bits;
  
%%
  %%%%%%%%%%%%%%%%%%%%%%% adding AWGN %%%%%%%%%%%%%%%%%%%%%%
  
  No = 100;
  deltas_r_time_awgn = zeros(size(deltas_r_time));
  AWGN = ((No./2)^.5) *randn(size(deltas_r_time));
  deltas_r_time_awgn = deltas_r_time + AWGN ;
  deltas_r_freq_awgn = fftshift(fft(deltas_r_time_awgn));
  %figure;
  %plot(t_axis,deltas_r_time_awgn);
  %title('after pulse shaper and awgn in time')
  %xlim([0 T_Bit*5])/
  %ylim([-1.5 1.5])
  %hold on;
  %figure;
  %plot(f_axis,deltas_r_freq_awgn);
  %title('after pulse shaper and awgn in frequency')
  %xlim([-100000 100000])
  %  ylim([0 10000])

%%
  %%%%%%%%%%%%%%%%%%%%%%% channel after AWGN %%%%%%%%%%%%%%%%%%%%%%
  
  deltas_out_channel_awgn_freq = channel .* deltas_r_freq_awgn;
  deltas_out_channel_awgn_time = ifft(ifftshift(deltas_out_channel_awgn_freq));
  %figure;
  %plot(t_axis,deltas_out_channel_awgn_time);
  %title('out of deltas + AWGN channel in time')
  %xlim([0  5*T_Bit])
  %ylim([-1.5 1.5])
  %figure;
  %plot(f_axis,deltas_out_channel_awgn_freq);
  %title('out of deltas + AWGN channel in freq')
  %xlim([-100000 100000])
  
%%
  %%%%%%%%%%%%%%%%%%%%%%% threshold with AWGN %%%%%%%%%%%%%%%%%%%%%%
  
  exp_bit_AWGN = zeros(1,N);
  for i=1:N_bits
      if deltas_out_channel_awgn_time(1,((i-1)*N_sq)+1) > 0+epsilon
          exp_bit_AWGN(1,i) = 1;
      else
          exp_bit_AWGN(1,i) = 0;
      end
  end
%%
  %%%%%%%%%%%%%%%%%%%%%%% compute BER with AWGN %%%%%%%%%%%%%%%%%%%%%%
  
  e_AWGN = 0;
  for i=1:N_bits
      if exp_bit_AWGN(1,i) ~= random_Bits(1,i)
          e_AWGN = e_AWGN+1;
      end
  end
  BER = e_AWGN/N_bits;
  disp(BER);
%%
%% part-2 Inter-Symbol Interference due to multi-path channels
%%
%%%%%%%%%%%% generation of Transmitted symbols (BPSK) %%%%%%%%%%%%%
L = 4000;     %No of paths
X_Tx = [];  
X_Tx = randi([0 1],[L 1]);  
A=length(X_Tx);
for i= 1 : A 
    if X_Tx(i) == 0
        X_Tx(i)= -1;
    end
end
%%
%%%%%%%%%%%% generation Matrix of channel coeffients %%%%%%%%%%%%%
coeff=MultipathChannel(L,1);
H= zeros(L,L);
for i =1:L
    for j = i:L
        H(j,i) = coeff(j - i +1);
    end 
end 
V = H * X_Tx ;
Eb_No_db = 0;       % The specified Eb/No value in dB
Energy_per_bit=1;
No=Energy_per_bit/( 10^(Eb_No_db/10) );
noise= randn(size(V))*sqrt(No/2); %generate Noise
Y = (H * X_Tx) + noise ; %getting the received signal Y
%%
%%%%%%%%%%%% Estimation of transmitted signal from Received signal %%%%%%%%%%%%%
Z= inv(H); %Equalize channel effect
X= Z * Y;  
X_Estimated =[]; 
%Decision Maker 
for i=1:A
    if X(i) > 0
        B =1; 
    else
        B=-1;
    end
    X_Estimated =[X_Estimated ; B];   %Estimated_Transmitted_signal
end
%Calculation of BER between transmitted symbols & estimated signal when
%variance = 1
BER = ComputeBER(X_Tx, X_Estimated);
%%
%%%%%%%%%%%% Estimation of BER vs Eb/No %%%%%%%%%%%%%
Eb_No_dB_vector = -15:0;
BER1=zeros(size(Eb_No_dB_vector));

for i= 1:length(Eb_No_dB_vector)
    
No=Energy_per_bit/( 10^(Eb_No_dB_vector(i)/10) );
noise= randn(size(V))*sqrt(No/2);
Y = (H * X_Tx) + noise ;
Z= inv(H);
X= Z * Y;
X_Estimated =[];
for j=1:A
    if X(j) > 0
        B =1; 
    else
        B =-1;
    end
X_Estimated =[X_Estimated ; B];
end
BER1(i) = ComputeBER(X_Tx, X_Estimated);
end
%Plotting BER vs Eb/No
figure
semilogy(Eb_No_dB_vector,BER1,'-xk','linewidth',2)
xlabel('Eb/No','linewidth',2)
ylabel('BER','linewidth',2)
%% part-3 Comparisons of coding techniques
%% Repetition code
%%
%%%%%%%%%%%Generates a sequence of bits %%%%%%%%%%%%%
p=.5;  % Channel parameter (probability of bit flipping)
n_bits=1000;         
bit_seq=randi([0,1],1,n_bits);
c_bits=5000;
r=n_bits/c_bits;
L=1./r;
%%
%%%%%%%%%%%%%%%%%%Generates a sequence of samples%%%%%%%%%%%
sample_seq=[];
for n=1:length(bit_seq)
c=repmat(bit_seq(n),[1,L]);
sample_seq=[sample_seq c];
end
%%
%%%%%%%%%BSC%%%%%%%%%%%%%%%%%
channel_effect = rand(size(sample_seq))<=p;
rec_sample_seq = xor(sample_seq,channel_effect);

%%
rec_bit_seq =[];
for n=1:L:length(rec_sample_seq)
    if (sum(rec_sample_seq(n:n+(L-1)))>(L/2))
    rec_bit_seq=[rec_bit_seq 1];
    else
    rec_bit_seq=[rec_bit_seq 0];
    end
end
%%
%%%%%%%%%%%% Compute  BER %%%%%%%%%%%%%%%
 a= (bit_seq == rec_bit_seq);
      E=0;
      for i=1:length(bit_seq)
         if a(i)==0 
             E=E+1;
         end
      end
BER=E/length(bit_seq);

%%
%%%%%%%%%%%% Effect of bit flipping probability on BER%%
p_vect          = 0:0.1:.5;              
BER_case_3_vec  = zeros(size(p_vect));  

for p_indx = 1:length(p_vect)
    
    channel_effect = rand(size(sample_seq))<=p_vect(p_indx);
    rec_sample_seq = xor(sample_seq,channel_effect);
   
    rec_bit_seq =[];
    for n=1:L:length(rec_sample_seq)
        if (sum(rec_sample_seq(n:n+(L-1)))>(L/2))
           rec_bit_seq=[rec_bit_seq 1];
        else
           rec_bit_seq=[rec_bit_seq 0];
        end
    end
    
     a= (bit_seq == rec_bit_seq);
      E=0;
      for i=1:length(bit_seq)
         if a(i)==0 
             E=E+1;
         end
      end
      
    BER_case_3_vec(p_indx)=E/length(bit_seq);
 
end
%%%


figure
plot(p_vect,BER_case_3_vec,'d-b','linewidth',2); hold on;
xlabel('Values of p','fontsize',10)
ylabel('BER','fontsize',10)
%% Convolutional code
%%

%% Polar codes
%%
%% Reliability sequence (Frozen positions)
Frozen_pos =[0 1 2 4 8 16 32 3 5 64 9 6 17 10 18 128 12 33 65 20 256 34 24 36 7 129 66 512 11 40 68 130 ...
   19 13 48 14 72 257 21 132 35 258 26 513 80 37 25 22 136 260 264 38 514 96 67 41 144 28 69 42 ...
   516 49 74 272 160 520 288 528 192 544 70 44 131 81 50 73 15 320 133 52 23 134 384 76 137 82 56 27 ...
   97 39 259 84 138 145 261 29 43 98 515 88 140 30 146 71 262 265 161 576 45 100 640 51 148 46 75 266 273 517 104 162 ...
   53 193 152 77 164 768 268 274 518 54 83 57 521 112 135 78 289 194 85 276 522 58 168 139 99 86 60 280 89 290 529 524 ...
   196 141 101 147 176 142 530 321 31 200 90 545 292 322 532 263 149 102 105 304 296 163 92 47 267 385 546 324 208 386 150 153 ...
   165 106 55 328 536 577 548 113 154 79 269 108 578 224 166 519 552 195 270 641 523 275 580 291 59 169 560 114 277 156 87 197 ...
   116 170 61 531 525 642 281 278 526 177 293 388 91 584 769 198 172 120 201 336 62 282 143 103 178 294 93 644 202 592 323 392 ...
   297 770 107 180 151 209 284 648 94 204 298 400 608 352 325 533 155 210 305 547 300 109 184 534 537 115 167 225 326 306 772 157 ...
   656 329 110 117 212 171 776 330 226 549 538 387 308 216 416 271 279 158 337 550 672 118 332 579 540 389 173 121 553 199 784 179 ...
   228 338 312 704 390 174 554 581 393 283 122 448 353 561 203 63 340 394 527 582 556 181 295 285 232 124 205 182 643 562 286 585 ...
   299 354 211 401 185 396 344 586 645 593 535 240 206 95 327 564 800 402 356 307 301 417 213 568 832 588 186 646 404 227 896 594 ...
   418 302 649 771 360 539 111 331 214 309 188 449 217 408 609 596 551 650 229 159 420 310 541 773 610 657 333 119 600 339 218 368 ...
   652 230 391 313 450 542 334 233 555 774 175 123 658 612 341 777 220 314 424 395 673 583 355 287 183 234 125 557 660 616 342 316 ...
   241 778 563 345 452 397 403 207 674 558 785 432 357 187 236 664 624 587 780 705 126 242 565 398 346 456 358 405 303 569 244 595 ...
   189 566 676 361 706 589 215 786 647 348 419 406 464 680 801 362 590 409 570 788 597 572 219 311 708 598 601 651 421 792 802 611 ...
   602 410 231 688 653 248 369 190 364 654 659 335 480 315 221 370 613 422 425 451 614 543 235 412 343 372 775 317 222 426 453 237 ...
   559 833 804 712 834 661 808 779 617 604 433 720 816 836 347 897 243 662 454 318 675 618 898 781 376 428 665 736 567 840 625 238 ...
   359 457 399 787 591 678 434 677 349 245 458 666 620 363 127 191 782 407 436 626 571 465 681 246 707 350 599 668 790 460 249 682 ...
   573 411 803 789 709 365 440 628 689 374 423 466 793 250 371 481 574 413 603 366 468 655 900 805 615 684 710 429 794 252 373 605 ...
   848 690 713 632 482 806 427 904 414 223 663 692 835 619 472 455 796 809 714 721 837 716 864 810 606 912 722 696 377 435 817 319 ...
   621 812 484 430 838 667 488 239 378 459 622 627 437 380 818 461 496 669 679 724 841 629 351 467 438 737 251 462 442 441 469 247 ...
   683 842 738 899 670 783 849 820 728 928 791 367 901 630 685 844 633 711 253 691 824 902 686 740 850 375 444 470 483 415 485 905 ...
   795 473 634 744 852 960 865 693 797 906 715 807 474 636 694 254 717 575 913 798 811 379 697 431 607 489 866 723 486 908 718 813 ...
   476 856 839 725 698 914 752 868 819 814 439 929 490 623 671 739 916 463 843 381 497 930 821 726 961 872 492 631 729 700 443 741 ...
   845 920 382 822 851 730 498 880 742 445 471 635 932 687 903 825 500 846 745 826 732 446 962 936 475 853 867 637 907 487 695 746 ...
   828 753 854 857 504 799 255 964 909 719 477 915 638 748 944 869 491 699 754 858 478 968 383 910 815 976 870 917 727 493 873 701 ...
   931 756 860 499 731 823 922 874 918 502 933 743 760 881 494 702 921 501 876 847 992 447 733 827 934 882 937 963 747 505 855 924 ...
   734 829 965 938 884 506 749 945 966 755 859 940 830 911 871 639 888 479 946 750 969 508 861 757 970 919 875 862 758 948 977 923 ...
   972 761 877 952 495 703 935 978 883 762 503 925 878 735 993 885 939 994 980 926 764 941 967 886 831 947 507 889 984 751 942 996 ...
   971 890 509 949 973 1000 892 950 863 759 1008 510 979 953 763 974 954 879 981 982 927 995 765 956 887 985 997 986 943 891 998 766 ...
   511 988 1001 951 1002 893 975 894 1009 955 1004 1010 957 983 958 987 1012 999 1016 767 989 1003 990 1005 959 1011 1013 895 1006 1014 1017 1018 ...
   991 1020 1007 1015 1019 1021 1022 1023]+1;

%% Polar Transformation Encoder
%%
%%%%%%%%%%%%Generates a sequence of bits %%%%%%%%%%%%%
n=4;
N=2^n;    
K=10;
bit_seq=randi([0,1],1,K);
%%
%%%%%%%%%%%%%%%%Generates Input signal %%%%%%%%%%%%%%%
R = Frozen_pos(Frozen_pos<=N); %reliability sequence for N
U=zeros(1,N);
U(R(N-K+1 : end))=bit_seq;
%%
%%%%%%%%%%Generates a Generator Matrix N*N %%%%%%%%%%%
%N=2^M
G2=[1 0 ; 1 1];
G=G2;
for j=1:(n-1) 
   G=kron(G2,G);  %Polar transformation
end
%%
%%%%%%%Encode Input Signal using Polar Transform%%%%%%%
X= U*G;
X_TX=[];
for i=1:length(X)
    if mod(X(i),2) == 0
     X_TX=[X_TX 0];
    else
     X_TX=[X_TX 1];   
    end
end
%% Channel Polarization
%%
%%%%%%%%%%%%BSC channel effect  %%%%%%%%%%%%%
p=.5;  % Channel parameter (probability of bit flipping)
channel_effect = rand(size(X_TX))<=p;
X_RX = xor(X_TX,channel_effect);



