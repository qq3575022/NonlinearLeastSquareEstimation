function J = getJacoN(y,x,N,class)

dx = 1e-12;
I = eye(length(x));
J = [ ];

if class == 1
    for j = 1 : length(x)
    e = I(:,j);
    E = getEP(y,x,N);
    Edx = getEP(y,x+dx*e,N);
    J(:,j) = (Edx-E)/dx;
    end
elseif class == 2
    for j = 1 : length(x)
    e = I(:,j);
    E = getEPV(y,x,N);
    Edx = getEPV(y,x+dx*e,N);
    J(:,j) = (Edx-E)/dx;
    end
elseif class == 3
    for j = 1 : length(x)
    e = I(:,j);
    E = getEPA(y,x,N);
    Edx = getEPA(y,x+dx*e,N);
    J(:,j) = (Edx-E)/dx;
    end
elseif class == 4
    for j = 1 : length(x)
    e = I(:,j);
    E = getEPVA(y,x,N);
    Edx = getEPVA(y,x+dx*e,N);
    J(:,j) = (Edx-E)/dx;
    end
elseif class == 5
    for j = 1 : length(x)
    e = I(:,j);
    E = getEPVAO(y,x,N);
    Edx = getEPVAO(y,x+dx*e,N);
    J(:,j) = (Edx-E)/dx;
    end
else
    for j = 1 : length(x)
    e = I(:,j);
    E = getEPAO(y,x,N);
    Edx = getEPAO(y,x+dx*e,N);
    J(:,j) = (Edx-E)/dx;
    end
end