clear ;
clc;
close all;
f1=[];
f2=[];
f3=[];
f4=[];
f5=[];
U=8;
alpha=zeros(448);
for k=1:448
    alpha(k , k)= exp(sqrt(-1).*2*pi*(k-1)./448);
end
l=-22.3:0.1:22.4;
a=0.4;
g1 = (sinc(l).*cos(pi*a*l))./((1-(2*a*l).^2)); %RC filter
n=0:1:447; %gfdm symbols
A=[]; %modulation matrix 
for i=1:448
    A=[A (circshift(g1,64.*fix((i-1)/64)).*exp(sqrt(-1)*2*pi*(fix((i-1)/64))))'];
end
A=A';
for z=1:1000
    %GFDM PAPR
    b= randi([0 1],1792,1);%binary vector
    y=[1:448]; 
    for i=1:1:448
        j=4*i-3;
        y(i)=b(j)*8+b(j+1)*4+b(j+2)*2+b(j+3); %decimal conversion
    end
    d=qammod(y,16);%data set vector
    d=d';
    x=A*d;
    pp=0;
    sp=0;
    for o=1:448 
        e=x(o).*conj(x(o));
        pp=max(pp,e);
        sp=sp+e;
    end
    pap=448*pp./sp;
    f1=[f1 10*log10(pap)];
    Beta=[];
    I=eye(448);
    PP_m=50;
    for j=1:8
        I=I*alpha;
        v=A*I*d;
        Beta=[Beta v];
    end
    gamma=Beta;
    for j=1:8
        for t=j+1:8
            h1=(1/sqrt(2)).*(Beta(1:448,j)+(sqrt(-1).*Beta(1:448,t)));
            h2=(1/sqrt(2)).*(Beta(1:448,j)-(sqrt(-1).*Beta(1:448,t)));
            gamma =[gamma h1 h2];
        end
    end
    for j=1:8
        pp=0;
        sp=0;
        v=Beta(1:448,j);
        v=v';
        for o=1:448
            e=v(o).*conj(v(o));
            pp=max(pp,e);
            sp=sp+e;
        end
        pap=448*pp./sp;
        pap =10*log10(pap);
        if(pap<PP_m)
            V=v;
            PP_m=pap;
        end
    end
    f2 =[f2 PP_m];
    PP_m=50;
    for j=1:64
        pp=0;
        sp=0;
        v=gamma(1:448,j);
        v=v';
        for o=1:448
            e=v(o).*conj(v(o));
            pp=max(pp,e);
            sp=sp+e;
        end
        pap=448*pp./sp;
        pap =10*log10(pap);
        if(pap<PP_m)
            V=v;
            PP_m=pap;
        end
        if(j==16)
            f3=[f3 PP_m];
        end
        if(j==32)
            f4=[f4 PP_m];
        end
    end
    f5=[f5 PP_m];
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
grid on;
xlabel('PAPR in dB');
ylabel('CCDF');
title('PAPR Reduction using proposed SLM technique');
legend('without slm','with slm 8U','with  proposed slm 16U ','with  proposed slm 32U','with  proposed slm 64U ');