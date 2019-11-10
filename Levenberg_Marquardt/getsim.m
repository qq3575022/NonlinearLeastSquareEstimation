function [phi_mu,phi_mod,r,rdot,diff,l] = getsim(x,f,Gt,M,X,PT,GT,GR,R,sigma,index,k,z,z_prev,T,l)

lambda = 3*10^8/f;

phi_prev_mod = mod(4*pi*sqrt((z_prev(1,k)-x(1))^2+(z_prev(2,k)-x(2))^2)/lambda,2*pi);
phi_mod      = mod(4*pi*sqrt((z(1,k)     -x(1))^2+(z(2,k)     -x(2))^2)/lambda,2*pi);
diff = phi_mod - phi_prev_mod;

if index == 1
    if (phi_prev_mod - phi_mod)> 2.1
        l = l + 1;
        phi_prev_mu = phi_prev_mod + (l-1)*2*pi;
        phi_mu = phi_mod + l*2*pi;

    elseif (phi_mod - phi_prev_mod)> 4.2
        l = l - 1;
        phi_prev_mu = phi_prev_mod + l*2*pi;
        phi_mu = phi_mod + (l-1)*2*pi;

    else
        phi_prev_mu = phi_prev_mod + l*2*pi;
        phi_mu = phi_mod + l*2*pi;
    end
    
elseif index == 2
    if k < 0.47*length(z)
        if (phi_prev_mod - phi_mod)> 2.9
            l = l + 1;
            phi_prev_mu = phi_prev_mod + (l-1)*2*pi;
            phi_mu = phi_mod + l*2*pi;

        else
            phi_prev_mu = phi_prev_mod + l*2*pi;
            phi_mu = phi_mod + l*2*pi;
        end
    else
        
        if (phi_mod - phi_prev_mod)> 3.0
            l = l - 1;
            phi_prev_mu = phi_prev_mod + l*2*pi;
            phi_mu = phi_mod + (l-1)*2*pi;

        else
            phi_prev_mu = phi_prev_mod + l*2*pi;
            phi_mu = phi_mod + l*2*pi;
        end
    end
    
else 
    if k < 0.5*length(z)
        if (phi_mod - phi_prev_mod)> 2.3
            l = l - 1;
            phi_prev_mu = phi_prev_mod + l*2*pi;
            phi_mu = phi_mod + (l-1)*2*pi;
            
        elseif (phi_prev_mod - phi_mod)> 5.0
            l = l + 1;
            phi_prev_mu = phi_prev_mod + (l-1)*2*pi;
            phi_mu = phi_mod + l*2*pi;
        else
            phi_prev_mu = phi_prev_mod + l*2*pi;
            phi_mu = phi_mod + l*2*pi;
        end
        
    else
        if (phi_prev_mod - phi_mod)> 3.0
            l = l + 1;
            phi_prev_mu = phi_prev_mod + (l-1)*2*pi;
            phi_mu = phi_mod + l*2*pi;
   
        else
            phi_prev_mu = phi_prev_mod + l*2*pi;
            phi_mu = phi_mod + l*2*pi;
        end
               
    end
end

phi_prev_conc = phi_prev_mu;  %icdf('Normal',uniphase,phi_prev_mu,sigma);
phi_conc      = phi_mu;       %icdf('Normal',uniphase2,phi_mu,sigma);

uniphase = rand;
delta_phi = phi_conc - phi_prev_conc;
delta_phi = icdf('Normal',uniphase,delta_phi,sigma);

H = sqrt((2*R*PT*GT*GR*Gt^2*lambda^4*X^2*M)/(4*pi*sqrt((z(1,k)-x(1))^2+(z(2,k)-x(2))^2))^4);
v = H;
unirand = rand;
% % if sigma>>v
% coeif = raylinv(unirand,sigma)
% v = raylinv(unirand,sigma); %

v = icdf('Normal',unirand,v,sigma); %sigma*randn + mu;

r   = ((2*PT*GT*GR*Gt^2*lambda^4*X^2*M*R)/((4*pi)^4*v^2))^(1/4);

rdot = lambda/(4*pi)*delta_phi*1/T;

end