b= randi([0 1],1024,1); %binary vector
y=[1:256]; 
for i=1:1:256
    j=4*i-3;
    y(i)=b(j)*8+b(j+1)*4+b(j+2)*2+b(j+3); %decimal conversion
end
d=qammod(y,16); %data set vector
l=-12.8:0.1:12.7;
%subplot(4,3,6);
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
[n, i]= meshgrid(n,i);
mesh(n,i,abs(A));
xlabel('sample index n');
ylabel('column index i');
zlabel('|A|');
title('GFDM Modulator Matrix');
xlim([0 260]);
ylim([0 260]);
zlim([0 1]);

