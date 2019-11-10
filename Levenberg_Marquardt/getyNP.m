function y = getyNP(r1,r2,r3,k,N)

y = NaN(3*N,1);

for i = 1:1:N
y(1+3*(i-1):3*i) = [r1(k+i-1);r2(k+i-1);r3(k+i-1)];
end

end