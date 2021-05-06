n=-5:0.001:5;
a=0.1;
g1 = (sinc(n).*cos(pi*a*n))./((1-(2*a*n).^2));
plot(n,g1);
hold on 
a=0.5;
g1 = (sinc(n).*cos(pi*a*n))./((1-(2*a*n).^2));
plot(n,g1);
hold on 
a=0.9;
g1 = (sinc(n).*cos(pi*a*n))./((1-(2*a*n).^2));
plot(n,g1);
xlabel('n');
ylabel('h(n)');
title('Impulse Response of RC filter with various roll-off factors');
legend('β=0.1','β=0.5','β=0.9');

