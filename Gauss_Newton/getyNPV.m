function y = getyNPV(r1,r1dot,r2,r2dot,r3,r3dot,k,N)

y = NaN(6*N,1);

for i = 1:1:N
y(1+6*(i-1):6*i) = [r1(k+i-1);r1dot(k+i-1);r2(k+i-1);r2dot(k+i-1);r3(k+i-1);r3dot(k+i-1)];
end

% y(1:6)   = [r1(k);  r1dot(k);  r2(k);  r2dot(k);  r3(k);  r3dot(k)];
% y(7:12)  = [r1(k+1);r1dot(k+1);r2(k+1);r2dot(k+1);r3(k+1);r3dot(k+1)];
% y(13:18) = [r1(k+2);r1dot(k+2);r2(k+2);r2dot(k+2);r3(k+2);r3dot(k+2)];
% y(19:24) = [r1(k+3);r1dot(k+3);r2(k+3);r2dot(k+3);r3(k+3);r3dot(k+3)];
% y(25:30) = [r1(k+4);r1dot(k+4);r2(k+4);r2dot(k+4);r3(k+4);r3dot(k+4)];
% y(31:36) = [r1(k+5);r1dot(k+5);r2(k+5);r2dot(k+5);r3(k+5);r3dot(k+5)];
% y(37:42) = [r1(k+6);r1dot(k+6);r2(k+6);r2dot(k+6);r3(k+6);r3dot(k+6)];
% y(43:48) = [r1(k+7);r1dot(k+7);r2(k+7);r2dot(k+7);r3(k+7);r3dot(k+7)];
% y(49:54) = [r1(k+8);r1dot(k+8);r2(k+8);r2dot(k+8);r3(k+8);r3dot(k+8)];
% y(55:60) = [r1(k+9);r1dot(k+9);r2(k+9);r2dot(k+9);r3(k+9);r3dot(k+9)];

end