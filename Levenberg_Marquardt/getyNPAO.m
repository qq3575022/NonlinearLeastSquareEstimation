function y = getyNPAO(r1,r2,r3,phi1,ax1,ax2,k,N)

y = NaN(6*N,1);

for i = 1:1:N
y(1+6*(i-1):6*i) = [r1(k+i-1);r2(k+i-1);r3(k+i-1);phi1(k+i-1);ax1(k+i-1);ax2(k+i-1)];
end

end