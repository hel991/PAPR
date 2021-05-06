l=6:0.1:31.5;
f1=[];
f2=[];
f3=[];
a1=0.1;
a2=0.4;
a3=0.9;
g1 = (sinc(l).*cos(pi*a1*l))./((1-(2*a1*l).^2)); %RC filter
g2 = (sinc(l).*cos(pi*a2*l))./((1-(2*a2*l).^2));
g3 = (sinc(l).*cos(pi*a3*l))./((1-(2*a3*l).^2));
n=0:1:255; %gfdm symbols
A=[]; %modulation matrix 
A2=[];
A3=[];
for i=1:256
    A=[A (circshift(g1,64.*fix((i-1)/64)).*exp(sqrt(-1)*2*pi*(fix((i-1)/64))))'];
    A2=[A2 (circshift(g2,64.*fix((i-1)/64)).*exp(sqrt(-1)*2*pi*(fix((i-1)/64))))'];
    A3=[A3 (circshift(g3,64.*fix((i-1)/64)).*exp(sqrt(-1)*2*pi*(fix((i-1)/64))))'];
end
A=A';
A2=A2';
A3=A3';
for z=1:20000
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
    pp1=0;
    sp1=0;
    for o=1:256
        e=x(o).*conj(x(o));
        pp1=max(pp1,e);
        sp1=sp1+e;
    end
    pap1=256*pp1./sp1;
    f1=[f1 10*log10(pap1)];
    x2=A2*d;
    pp2=0;
    sp2=0;
    for o=1:256
        e2=x2(o).*conj(x2(o));
        pp2=max(pp2,e2);
        sp2=sp2+e2;
    end
    pap2=256*pp2./sp2;
    f2=[f2 10*log10(pap2)];
    x3=A3*d;
    pp3=0;
    sp3=0;
    for o=1:256
        e3=x3(o).*conj(x3(o));
        pp3=max(pp3,e3);
        sp3=sp3+e3;
    end
    pap3=256*pp3./sp3;
    f3=[f3 10*log10(pap3)];
end
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
xlim([0 7]);
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
xlim([0 7]);
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
xlim([0 7]);
ylim([10^-2 10^0]);
grid on;
xlabel('PAPR in dB');
ylabel('CCDF');
title('CCDF OF GFDM WITH DIFFERENT ROLL-OFF FACTORs');
legend('β=0.1','β=0.5','β=0.9');
