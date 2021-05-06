clear, clc, close all;
f=[];
h=[];
l=-12.8:0.1:12.7;
a=0.9;
g1 = (sinc(l).*cos(pi*a*l))./((1-(2*a*l).^2)); %RC filter
n=0:1:255; %gfdm symbols
A=[]; %modulation matrix 
for i=1:256
    A=[A (circshift(g1,64.*fix((i-1)/64)).*exp(sqrt(-1)*2*pi*(fix((i-1)/64))))'];
end
A=A';
i=0:1:255;
n=0:1:255;
k=1:1:64;
c1=[];
for t=-32:1:31
    c1=[c1 exp(sqrt(-1).*2*pi.*k.*t/64)'];
end
for z=1:1000
    %GFDM PAPR
    b= randi([0 1],1024,1);%binary vector
    y=[1:256]; 
    for i=1:1:256
        j=4*i-3;
        y(i)=b(j)*8+b(j+1)*4+b(j+2)*2+b(j+3); %decimal conversion
    end
    d=qammod(y,16); %data set vector
    x=(A)*d';
    o=1:1:256; 
    pp=0;
    sp=0;
    for o=1:256 
        e=x(o).*conj(x(o));
        pp=max(pp,e);
        sp=sp+e;
    end
    pap=256*pp./sp;
    f=[f 10*log10(pap)];
    %OFDM PAPR
    x21=d(1:64)*c1;
    x22=d(65:128)*c1;
    x23=d(129:192)*c1;
    x24=d(193:256)*c1;
    x2=[x21 x22 x23 x24];
    o1=1:1:256; 
    pp1=0;
    sp1=0;
    for o1=1:256
        e1=x2(o1).*conj(x2(o1));
        pp1=max(pp1,e1);
        sp1=sp1+e1;
    end
    pap1=256*pp1./sp1;
    h=[h 10*log10(pap1)];
end
%subplot(2,1,1);
%stem(f);
Th = 20*log10([1:0.1:10]);
kk = 1;
th_current = Th(1);
z = sort(f);
L = -inf*ones(1, length(Th));
for ii = 1:(length(z)-1)
 if (z(ii) > th_current)
 L(kk) = ii - 1;
 kk = kk + 1;
 th_current = Th(kk);
 end
end
CCDF = (length(z) - L)/length(z);
semilogy(Th, CCDF);
xlim([0 12]);
ylim([10^-2 10^0]);
grid on;
hold on;
%subplot(2,1,2);
%stem(h);
Th = 20*log10([1:0.1:10]);
kk = 1;
th_current = Th(1);
z = sort(h);
L = -inf*ones(1, length(Th));
for ii = 1:(length(z)-1)
 if (z(ii) > th_current)
 L(kk) = ii - 1;
 kk = kk + 1;
 th_current = Th(kk);
 end
end
CCDF = (length(z) - L)/length(z);
semilogy(Th, CCDF);
xlim([0 12]);
ylim([10^-2 10^0]);
grid on;
xlabel('PAPR in dB');
ylabel('CCDF ');
title('CCDF COMPARISON BETWEEN GFDM AND OFDM');
legend('GFDM','OFDM');