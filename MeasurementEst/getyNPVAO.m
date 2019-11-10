function y = getyNPVAO(r1,r1dot,r2,r2dot,r3,r3dot,phi1,ax1,ax2,k,N)

y = NaN(9*N,1);

for i = 1:1:N
y(1+9*(i-1):9*i) = [r1(k+i-1);r1dot(k+i-1);r2(k+i-1);r2dot(k+i-1);r3(k+i-1);r3dot(k+i-1);phi1(k+i-1);ax1(k+i-1);ax2(k+i-1)];
end

end