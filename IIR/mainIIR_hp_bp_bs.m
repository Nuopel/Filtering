% 
%  File    :   mainIIR_hp_bp_bs
% 
% Author1: Samuel Dupont
% Date:    November 2016 
%   
%  Description :  Explore IIR filtering for high pass departing from analitical
%                   formulation in order to produce the filter coefficients.
%                   It departs from the analog coefficient, then using the
%                   bilinear transform it produces the digital numbers;
%                   Contrary to MainIIR.m it use a straight Biquad
%                   implementation
%                 
%
%  Reference : 
%                 https://github.com/ruohoruotsi/Butterworth-Filter-Design/blob/master/Butterworth.cpp
%                 http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
%                 https://en.wikipedia.org/wiki/Butterworth_filter
%                 https://en.wikipedia.org/wiki/Bilinear_transform
    clear variables; close all; clc

    %% Define the various constant of the filter
    c.FS = 44100;
    c.TS = 1/c.FS;
    c.NY = c.FS/2;
    
    % Filt property
    c.FC = 19000;%c.NY/2;%cut off freq
    c.ORDER = 1;
    v.k = (1:c.ORDER).';% vector 1,2...order of the filter
    
    %warp filter
    c.WC= 2/c.TS * tan( pi * c.FC / c.FS);
    K=2/(c.TS*c.WC);
    K2=K^2;
    
    %% Create the prototype continuous filter
    Coef=zeros(round(c.ORDER/2),6);
    fin=1:round(c.ORDER/2)-1;
    last=round(c.ORDER/2);

    Coef(fin,1) = 1;
    Coef(fin,2) = 0;
    Coef(fin,3) = 0;
    Coef(fin,4) = 1;
    Coef(fin,5) = -2*cos(pi*(2*fin+c.ORDER-1)/2/c.ORDER);
    Coef(fin,6) = 1;

    if mod(c.ORDER,2)==0
        Coef(end,1) = 1;
        Coef(end,2) = 0;
        Coef(end,3) = 0;
        Coef(end,4) = 1;
        Coef(end,5) = -2*cos(pi*(2*last+c.ORDER-1)/2/c.ORDER);
        Coef(end,6) = 1;
    else
        Coef(end,1) = 1;
        Coef(end,2) = 0;
        Coef(end,3) = 0;
        Coef(end,4) = 1;
        Coef(end,5) = 1;
        Coef(end,6) = 0;
    end

%% Calculate digital filter coefficients
    Coefd=zeros(last,6);
    
    Coefd(fin,1) = 1;
    Coefd(fin,2) = -2*K2./ ( K2 );
    Coefd(fin,3) = 1;
    Coefd(fin,4) = 1;
    Coefd(fin,5) = ( 2 - 2*K2 ) ./ ( K2 + Coef(fin,5)*K + 1 );
    Coefd(fin,6) = ( K2 -Coef(fin,5)*K + 1 ) ./ ( K2 + Coef(fin,5)*K + 1 );
    
    Gz = prod(K2./( K2 + Coef(fin,5)*K + 1));
    
    if mod(c.ORDER,2)==1
        Coefd(last,1) = Gz*K / ( Coef(last,5) +K );
        Coefd(last,2) = -Gz*K / ( Coef(last,5) +K );
        Coefd(last,4) = 1;
        Coefd(last,5) = ( +Coef(last,5) -K ) / ( Coef(last,5) +K );
    else
        Coefd(last,1) = Gz*K2./ ( K2 + Coef(last,5)*K + 1 );
        Coefd(last,2) = -Gz*2*K2./ ( K2 + Coef(last,5)*K + 1 );
        Coefd(last,3) = Gz*K2./ ( K2 + Coef(last,5)*K + 1 );
        Coefd(last,4) = 1;
        Coefd(last,5) = ( 2 - 2*K2 ) ./ ( K2 + Coef(last,5)*K + 1 );
        Coefd(last,6) = ( K2 -Coef(last,5)*K + 1 ) ./ ( K2 + Coef(last,5)*K + 1 );
    end

%% verif
    freqz(Coefd(1:3),Coefd(4:6),10000);
    hold on
    [b,a]=butter(c.ORDER,c.FC/c.NY,'high');sos=tf2sos(b,a);freqz(b,a)