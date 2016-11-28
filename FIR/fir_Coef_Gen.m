% 
%  File    :   fir_Coef_Gen  
% 
%  Author1 :   Samuel Dupont
%  Date    :   November 2016
% 
%  Course  :   Internet peregrination
%   
%  Description :   
%                     
%    1) This program shows how to produce FIR filter coefficients using sinus cardinal signal 
%       It emcompass also the use of some windows functions (Black Man) to reduce the ripple 
%       of such method.      
%
%    2) Data-structure : array, struct.
%                     
%    3) What algorithms, techniques, etc, are used in this program 
% 
%  Generate plot of the  high pass, low pass and merge of the both filter
%  You can choose the parameters of the filter in the constant part
%  definition (sampling freq, number of points in the filters, cut off frequency...)
%                     
% 
clear variables; close all;clc
tic
%% constants part definition
    ct.FS = 44100;
    ct.NY = ct.FS/2;
    ct.NFFT = 2048;
    
    fir.NT = 31;%number of points of the filters
    fir.BW = 4410/ct.NY;%cut off freq of the filter
    fir.BW1 =   0.3;% cut off freq 1 bandpass filter
    fir.BW2 = 0.6;%cut off freq 2 bandpass filter
    

    
%% sinc fir implementation
    if mod(fir.NT/2,1) == 0
        fir.t = 1-fir.NT/2 : fir.NT/2;%index of the filter 
        fir.n = 1:fir.NT; %Rect index of abscisse, handle odd or even fir.NT 
    else
        fir.t = ceil(-fir.NT/2+0.5) : floor(fir.NT/2);
        fir.n = 0.5:fir.NT-0.5; %Rect index of abscisse, handle odd or even fir.NT 
    end
    fir.freq=FreqVect(ct.FS,ct.NFFT);
    
    fir.BlackMan = 0.42 - 0.5*cos(2*pi*fir.n/fir.NT) + 0.08*cos(4*pi*fir.n/fir.NT) ;% Blackman windows 
    
    fir.lowPass = sinc(pi.*fir.t*fir.BW/pi).*fir.BW ;% generate filter
    fir.lowPassBM = fir.BlackMan .* fir.lowPass; %windowed sinc
 
    fir.bandpass = sinc(pi*fir.BW2*(fir.t)/pi) .* fir.BW2 ...
                   - sinc(pi*fir.BW1*(fir.t)/pi) .* fir.BW1 ;
    fir.bandpassBM =  fir.bandpass .* fir.BlackMan;  

    fir.highPass = -fir.lowPass;
    fir.highPass(round(fir.NT/2))=1+fir.highPass(round(fir.NT/2));%% spectral inversion , handle odd or even NT
    fir.highPassBM = fir.highPass .* fir.BlackMan;  
toc
 %%  frequency shift : mirroring by fs/2
%      test=ones(1,fir.NT);
%      test(2:2:end)=-test(2:2:end);
%      fir.lowPass = fir.lowPassBM.*test;%--->shift low pass to highPass  by half fs,
%                                         %%cut off being fs/2-bw
%      fir.highPass = fir.highPassBM.*test;
                                        
%      fir.lowPassBM = fir.lowPassBM.*test;
 %%  frequency shift : move the filter by X Hz       
%     mv=8000/ct.NY;% move 4000 Hz
%     test= cos(fir.n*pi*mv); %filter= lowpass shifted by 1000Hz
%      fir.lowPass = fir.lowPass.*test;%--->shift low pass to highPass  by half fs,
%                                         %cut off being fs/2-bw
%      fir.lowPassBM = fir.lowPassBM.*test;
                                        
%% plot   
tic
    %low pass
    figure(1)
        subplot(211);
            plot(fir.t,fir.lowPass);
               % title('Low pass sinc filters')
                xlabel('Sample');ylabel('Amplitude')
            hold on;
            plot(fir.t,fir.lowPassBM);
            grid on ;
                legend('Rect','Blackman')
        subplot(212)
            plot(fir.freq,db(abs(fft(fir.lowPass,ct.NFFT)))); 
            hold on;         
            plot(fir.freq,db(abs(fft(fir.lowPassBM,ct.NFFT))));
            plot([fir.BW*ct.NY fir.BW*ct.NY],[min(db(abs(fft(fir.lowPassBM,ct.NFFT))))...
                     max(db(abs(fft(fir.lowPassBM,ct.NFFT))))],'g');  
            
                xlabel('Frequency [Hz]');ylabel('Amplitude [dB]')
                xlim([0 ct.NY]);    
            grid on
                legend('Rect','Blackman')
    
    %high pass
    figure(2)
        subplot(211);
                title('High pass sinc filters')
            plot(fir.t,fir.highPass);
                xlabel('Sample');ylabel('Amplitude')
                xlim([fir.t(1) fir.t(end)]);
            hold on;
            plot(fir.t,fir.highPassBM);
            grid on ;
                legend('Rect','Blackman')
        subplot(212)
            plot(fir.freq,db(abs(fft(fir.highPass,ct.NFFT))));
            hold on;
            plot(fir.freq,db(abs(fft(fir.highPassBM,ct.NFFT))));
                xlabel('Frequency [Hz]');ylabel('Amplitude[dB]')
                xlim([0 ct.NY]); 
                plot([fir.BW*ct.NY fir.BW*ct.NY],[min(db(abs(fft(fir.lowPassBM,ct.NFFT))))...
                      max(db(abs(fft(fir.lowPassBM,ct.NFFT))))],'g');   
            grid on
                legend('Rect','Blackman')
 
                    %high pass
    figure(3)
        subplot(211);
                title('Bandpass sinc filters')
            plot(fir.t,fir.bandpass);
                xlabel('Sample');ylabel('Amplitude')
                xlim([fir.t(1) fir.t(end)]);
            hold on;
            plot(fir.t,fir.bandpassBM);
            grid on ;
                legend('Rect','Blackman')
        subplot(212)
            plot(fir.freq,db(abs(fft(fir.bandpass,ct.NFFT))));
            hold on;
            plot(fir.freq,db(abs(fft(fir.bandpassBM,ct.NFFT))));
                xlabel('Frequency [Hz]');ylabel('Amplitude[dB]')
                xlim([0 ct.NY]); 
            plot([fir.BW1*ct.NY fir.BW1*ct.NY],[min(db(abs(fft(fir.bandpassBM,ct.NFFT))))...
                      max(db(abs(fft(fir.bandpassBM,ct.NFFT))))],'g');   
            plot([fir.BW2*ct.NY fir.BW2*ct.NY],[min(db(abs(fft(fir.bandpassBM,ct.NFFT))))...
                      max(db(abs(fft(fir.bandpassBM,ct.NFFT))))],'g');       
            grid on
             legend('Rect','Blackman')
    %merge
    figure(4)
        semilogx(fir.freq,db(abs(fft(fir.lowPass,ct.NFFT))+abs(fft(fir.highPass,ct.NFFT))));
            xlabel('Sample');ylabel('Amplitude')
            xlim([fir.t(1) fir.t(end)]);
        hold on;
        semilogx(fir.freq,db(abs(fft(fir.lowPassBM,ct.NFFT))+abs(fft(fir.highPassBM,ct.NFFT))));
            xlabel('Frequency [Hz]');ylabel('Amplitude [dB]')
            xlim([0 ct.NY]);    
            %ylim([-1 max(db(abs(fft(fir.lowPass))+abs(fft(fir.highPass))))]);
        grid on
        legend('Rect','Blackman');
        title('Merging of sinc filters')

toc
%% windows
    win.rect=ones(fir.NT,1);
    win.bart=1-2*abs(fir.n.'-fir.NT/2)/fir.NT;
    win.blackman=0.42 - 0.5*cos(2*pi*fir.n.'/fir.NT) + 0.08*cos(4*pi*fir.n.'/fir.NT) ;% Blackman windows 
    win.hamming=0.54-0.46*cos(2*pi*fir.n.'/fir.NT);
    win.hanning=0.5-0.5*cos(2*pi*fir.n.'/fir.NT);
    
    figure(5)
        plot(fir.n,[win.rect, win.bart,win.blackman,win.hamming,win.hanning])
        ylim([0 1.1])
        xlim([fir.n(1) fir.n(end)])
        legend('Rect','Bartlett','Blackman','Hamming','Hanning')
    
    figure(6)
        plot(fir.freq,db(abs(fft(fir.lowPass.*win.rect.',ct.NFFT))));
            hold on 
        plot(fir.freq,db(abs(fft(fir.lowPass.*win.bart.',ct.NFFT))));
        plot(fir.freq,db(abs(fft(fir.lowPass.*win.blackman.',ct.NFFT))));
        plot(fir.freq,db(abs(fft(fir.lowPass.*win.hamming.',ct.NFFT))));
        plot(fir.freq,db(abs(fft(fir.lowPass.*win.hanning.',ct.NFFT))));
            legend('Rect','Bartlett','Blackman','Hamming','Hanning')
             xlabel('Frequency [Hz]');ylabel('Amplitude [dB]')
            ylim([-100 1])
            xlim([20 ct.NY])

%% point number influence
    fir.BW=0.5;
    fir.NT=[7; 30; 240];
    figure(7)
    for i=1:length(fir.NT)
        if mod(fir.NT(i)/2,1) == 0
            fir.t = 1-fir.NT(i)/2 : fir.NT(i)/2;%index of the filter 
            fir.n = 1:fir.NT(i); %Rect index of abscisse, handle odd or even fir.NT(i) 
        else
            fir.t = ceil(-fir.NT(i)/2+0.5) : floor(fir.NT(i)/2);
            fir.n = 0.5:fir.NT(i)-0.5; %Rect index of abscisse, handle odd or even fir.NT(i) 
        end
        length(fir.t)

        subplot(211)
        fir.lowPass = sinc(pi.*fir.t*fir.BW/pi).*fir.BW ;% generate filter
        plot(fir.freq,(abs(fft(fir.lowPass,ct.NFFT)))); 
            xlabel('Frequency [Hz]');ylabel('Amplitude [lin]')
            xlim([0 ct.NY]);      
            hold on;
            grid on;

        subplot(212)
        plot(fir.freq,db(abs(fft(fir.lowPass,ct.NFFT)))); 
            xlabel('Frequency [Hz]');ylabel('Amplitude [dB]')
            xlim([0 ct.NY]);  
            ylim([-60 10]);
            hold on;
            grid on;

    end 
    subplot(211)
    legend(num2str(fir.NT))
