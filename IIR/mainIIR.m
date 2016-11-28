% 
%  File    :   mainIIR
% 
% Author1: Samuel Dupont
% Date:    November 2016 
%   
%  Description :  Explore IIR filtering departing from analitical
%                   formulation in order to produce the filter coefficients.
%                   It departs from the analog coefficient, then using the
%                   bilinear transform it produces the digital numbers;
% 
%                   The value of the filter can be modified in the first
%                   part of the program.
%                 
%
%  Reference : 
%                 https://github.com/ruohoruotsi/Butterworth-Filter-Design/blob/master/Butterworth.cpp
%                 http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
%                 https://en.wikipedia.org/wiki/Butterworth_filter
%                 https://en.wikipedia.org/wiki/Bilinear_transform
    clear variables; close all; clc

%% Constant definition part
    c.FS = 44100;
    c.TS = 1/c.FS;
    c.NY = c.FS/2;
    
    % Filt property
    c.FC = 3000;%c.NY/2;%cut off freq
    c.WC = 2*pi*c.FC;% FC=freq*2/Fs
    c.ORDER = 4;
    v.k = (1:c.ORDER).';% vector 1,2...order of the filter
 %% Generate the butterworth coefficients in continous domain

    f.pk =-c.WC * exp(1i.* (2.*v.k+c.ORDER-1)*pi / (2*c.ORDER));% pole of the butterworth filt
    f.gs =c.WC^c.ORDER;% gain
    
    [b,a]=zp2tf([], f.pk, f.gs); % convert to transfert function
    [p.h,p.w]=freqs(b,a,2048); % continuous filter

    %plot
    figure(1)
    subplot(211)
        semilogx(p.w/(2*pi),db(abs(p.h)))
        grid on
        p.dyn=-min(db(abs(p.h)))+ round(max(db(abs(p.h))));
        axis([p.w(1)/(2*pi) p.w(end)/(2*pi) min(db(abs(p.h))) round(max(db(abs(p.h))))+0.05*p.dyn ])
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [dB]')
        hold on
        plot([c.FC c.FC],[min(db(abs(p.h))) round(max(db(abs(p.h))))],'g')
        plot([p.w(1)/(2*pi) p.w(end)/(2*pi) ],[-3 -3],'g')
        text(c.FC,-3,'\leftarrow Fc = -3dB')
     subplot(212)
        semilogx(p.w/(2*pi),unwrap(angle(p.h)))
        xlabel('Frequency [Hz]');
        ylabel('Angle [rad]');
        xlim([p.w(1)/(2*pi) p.w(end)/(2*pi)])
        grid on
%% Change to z domain using bilinear transform with warping method 1    
    c.W0= 2/c.TS * tan( pi * c.FC / c.FS);% warped cut off freq
    f.pk =-exp(1i.* (2.*v.k+c.ORDER-1)*pi / (2*c.ORDER));% normalised pole
    c.B=2./(c.TS*c.W0);% n front of bilinear transform (with warped coef)
    
    f.pkz2 = (c.B-f.pk) ./ (c.B+f.pk);% discrete pole
    f.gz2 = (1./prod(c.B+f.pk));% gain
    f.zk = ones(c.ORDER,1);% zeros

    [b,a] = zp2tf(f.zk,f.pkz2,f.gz2);%transform zeros poles form into transfert function (b,a)
    
    % plot
    figure(2)
    
    subplot(211)
    [p.h2,p.w2]=freqz(b,a);%obtain bode from b,a coef
        plot(p.w2/pi,db(abs(p.h2)))
        grid on
        p.dyn=-min(db(abs(p.h2)))+ round(max(db(abs(p.h2))));
        axis([p.w2(1)/pi p.w2(end)/pi min(db(abs(p.h2))) round(max(db(abs(p.h2))))+p.dyn*0.05]);
        xlabel('Frequency [Hz]');
        ylabel('Amplitude [dB]');
        hold on
        plot([c.FC/c.NY c.FC/c.NY],[min(db(abs(p.h2))) round(max(db(abs(p.h2))))+p.dyn*0.05],'g')
        text(c.FC/c.NY,-3,'\leftarrow Fc = -3dB')
          
    subplot(212)
        plot(p.w2/(2*pi),unwrap(angle(p.h2)))
        xlabel('Frequency [Hz]');
        ylabel('Angle [rad]');
        xlim([p.w2(1)/(2*pi) p.w2(end)/(2*pi)])
        grid on   
        
%% Difference between warped and normal
    opt=0;
    if opt==1
        c.FC=[3000  15000];
        c.WC=2*pi*c.FC;
        for ii=1:length(c.WC)
            % warped coefficient (like before)
            c.W0= 2/c.TS * tan( c.WC(ii)/2 / c.FS);
            f.pk =-exp(1i.* (2.*v.k+c.ORDER-1)*pi / (2*c.ORDER));
            c.B=2./(c.TS*c.W0);
            f.pkz2 = (c.B-f.pk) ./ (c.B+f.pk);
            f.gz2 = (1./prod(c.B+f.pk));

            [b,a] = zp2tf(f.zk,f.pkz2,f.gz2);
            [p.h2,p.w2]=freqz(b,a);

            % normal (without warped coefficient)
            f.pkz2 = (c.B-f.pk) ./ (c.B+f.pk);
            f.gz2 = (1./prod(c.B+f.pk));

            c.B=2./(c.TS*c.WC(ii));
            f.pkz3 = (c.B-f.pk) ./ (c.B+f.pk);
            f.gz3 = (1./prod(c.B+f.pk));

            [b,a] = zp2tf(f.zk,f.pkz3,f.gz3);
            [p.h3,p.w3]=freqz(b,a);

            %plot
            figure(4)

                plot(p.w2/pi,db(abs(p.h2)),'b')
                hold on
                grid on
                plot(p.w3/pi,db(abs(p.h3)),'r--')
            xlabel('Frequency [Hz]');
            ylabel('Amplitude [dB]');
            text(c.FC(ii)/c.NY,-3,'\leftarrow Fc = -3dB')
            plot([c.FC(ii)/c.NY c.FC(ii)/c.NY],[min(db(abs(p.h2))) round(max(db(abs(p.h2))))+p.dyn*0.05],'g')
            axis([p.w2(1)/pi p.w2(end)/pi min(db(abs(p.h2))) round(max(db(abs(p.h2))))+p.dyn*0.05]);
        end
           legend('warped', 'unwarped','location','southwest')
    end