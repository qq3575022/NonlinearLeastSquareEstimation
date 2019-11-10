function y = getyNPV(r1,r1dot,r2,r2dot,r3,r3dot,k,N)

y = NaN(6*N,1);

for i = 1:1:N
y(1+6*(i-1):6*i) = [r1(k+i-1);r1dot(k+i-1);r2(k+i-1);r2dot(k+i-1);r3(k+i-1);r3dot(k+i-1)];
end

end