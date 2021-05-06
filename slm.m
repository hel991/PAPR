clear ,clc, close all;
f1=[];
f2=[];
f3=[];
f4=[];
f5=[];
seq =zeros(256,64);
for R=1:64
    temp=[];
    temp = zadoffChuSeq(R,257);
    seq(1:256,R)=temp(2:257);
end
gamma=seq;
l=-12.8:0.1:12.7;
a=0.9;
g1 = (sinc(l).*cos(pi*a*l))./((1-(2*a*l).^2)); %RC filter
n=0:1:255; %gfdm symbols
A=[]; %modulation matrix 
for i=1:256
    A=[A (circshift(g1,64.*fix((i-1)/64)).*exp(sqrt(-1)*2*pi*(fix((i-1)/64))))'];
end
A=A';
for z=1:10000
    %GFDM PAPR
    b= randi([0 1],1024,1);%binary vector
    y=[1:256]; 
    for i=1:1:256
        j=4*i-3;
        y(i)=b(j)*8+b(j+1)*4+b(j+2)*2+b(j+3); %decimal conversion
    end
    d=qammod(y,16);%data set vector
    d=d';
    x=A*d;
    pp=0;
    sp=0;
    for o=1:256
        e=x(o).*conj(x(o));
        pp=max(pp,e);
        sp=sp+e;
    end
    pap=256*pp./sp;
    f1=[f1 10*log10(pap)];
    PP_m=50;
for j=1:64
    v=A*(gamma(1:256,j).*d);
    v=v';
    o=1:1:256; 
    pp=0;
    sp=0;
    for o=1:256 
        e=v(o).*conj(v(o));
        pp=max(pp,e);
        sp=sp+e;
    end
    pap=256*pp./sp;
    pap =10*log10(pap);
    if(pap<PP_m)
        V=v; 
        PP_m=pap;
    end
    if(j==8)
        f2=[f2 PP_m];
    end
    if(j==16)
        f3=[f3 PP_m];
    end
    if(j==32)
        f4=[f4 PP_m];
    end
end
f5 =[f5 PP_m];
end
% subplot(2,1,1);
% stem(f1);
% subplot(2,1,2);
% stem(f2);
Th = 20*log10([1:0.1:10]);
kk = 1;
th_current = Th(1);
z = sort(f1);
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
xlim([0 10]);
ylim([10^-2 10^0]);
grid on;
hold on;
Th = 20*log10([1:0.1:10]);
kk = 1;
th_current = Th(1);
z = sort(f2);
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
xlim([0 10]);
ylim([10^-2 10^0]);
grid on;
hold on;
Th = 20*log10([1:0.1:10]);
kk = 1;
th_current = Th(1);
z = sort(f3);
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
xlim([0 10]);
ylim([10^-2 10^0]);
grid on;
hold on;
Th = 20*log10([1:0.1:10]);
kk = 1;
th_current = Th(1);
z = sort(f4);
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
xlim([0 10]);
ylim([10^-2 10^0]);
grid on;
hold on;
Th = 20*log10([1:0.1:10]);
kk = 1;
th_current = Th(1);
z = sort(f5);
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
xlim([0 10]);
ylim([10^-2 10^0]);
grid on;
xlabel('PAPR in dB');
ylabel('CCDF');
title('PAPR Reduction using SLM(using Zhadoff chu phase seq) technique');
legend('original gfdm','using slm with U=8 ','using slm with U=16','using slm U=32','using slm with U=64 ');