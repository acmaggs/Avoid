clear
load jl
figure(2);
histogram(jl(:,1))
title('left')
figure(1);plot(jl(:,1), jl(:,2),'x')
title('left');
x=.8:.01:1.5;
y=1 ./ x.^12 - 1 ./ x.^6 +  1.58376 * x;
figure(3)
plot(x,y)


load jr
figure(20);
histogram(jr(:,1))
title('right')

figure(21);plot(jr(:,1), jr(:,2),'x')
title('right');
