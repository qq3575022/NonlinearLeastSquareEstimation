function E = getEP(y,x,N)

global T;
HN = NaN(3*N,1);

Gain = [1,T,0,0,0,0;
        0,1,T,0,0,0;
        0,0,1,0,0,0;
        0,0,0,1,T,0;
        0,0,0,0,1,T;
        0,0,0,0,0,1];

for i = 1:1:N
HN(1+3*(i-1):3*i)  = getHxkP(Gain^(i-1)*x);%C*getHxkPVA(Gain^(i-1)*x);
end

E = y - HN;

end