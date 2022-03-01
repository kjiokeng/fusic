%retourne R, theta0, theta1, t0 et t1 des l0 et l1 orth avec posMur=(a,0)
%pour x=a et posMur=(0,b) pour y=b
function [R, theta0,theta1,t0,t1] = getPaths(posC, posMur)
    if(posMur(1)==0)
        xr = posMur(2)*posC(1)/(2*posMur(2)-posC(2));
        yr = posMur(2); 
    else
        yr=posC(2)/( 1+ (posMur(1)-posC(1)) / posMur(1) );
        xr = posMur(1);
    end
    theta1= atand(xr/yr); %mettre moins avant d'entree cette valeur dans Demleloc, pareil pour en bas
    theta0= atand(posC(1)/posC(2));
    t1= norm([xr,yr])/(3*10^8);
    t0= norm(posC)/(3*10^8);
    R=[xr,yr];
end