function E = getEPAO(y,x,N)
global T

HN = NaN(6*N,1);

Gain = [1,T,0,0,0,0,0,0;
        0,1,T,0,0,0,0,0;
        0,0,1,0,0,0,0,0;
        0,0,0,1,T,0,0,0;
        0,0,0,0,1,T,0,0;
        0,0,0,0,0,1,0,0;
        0,0,0,0,0,0,1,T;
        0,0,0,0,0,0,0,1];

for i = 1:1:N

HN(1+6*(i-1):6*i)  = getHxkPAO(Gain^(i-1)*x);

end

E = y - HN;

end