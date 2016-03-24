function [ h, mag_h, i, omega, e_vector, mag_e, w, theta2, rp, zp, epslon ] = Gauss ( rho1_hat, rho2_hat, rho3_hat, R1, R2, R3, t1, t2, t3, muo, f0, R, improvement )
% this function is about Gauss method of preliminary orbit determination
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% inputs:
% rho1_hat               : direction cosine vectors of observer 1 (1x3)
% rho2_hat               : direction cosine vectors of observer 2 (1x3)
% rho3_hat               : direction cosine vectors of observer 3 (1x3)
% R1                         : position of observer 1 (1x3)
% R2                         : position of observer 2 (1x3)
% R3                         : position of observer 3 (1x3)
% t1                          : time of observer 1 in sec
% t2                          : time of observer 2 in sec
% t3                          : time of observer 3 in sec
% muo                      : Gravitational Parameter
% f0                          : initial value of r2
% R                           : earth raduis
% improvement       : if 0 no improve in orbital elements, if 1 improve in orbital elements
%% outputs:
% h              : specific angular momentum vector
% mag_h     : specific angular momentum magnitude
% i               : inclination angle in degree
% omega     : right ascension of the ascending node in degree
% e              : eccentricity vector
% mag_e     : eccentricity magnitude
% w             : argument of perigee in degree
% theta2       : true anomaly of r2 in degree
% epslon      :  specific energy
% ---------------------------------------------------------------------------------------------------------------------------------------------------------
tau1=t1-t2;
tau3=t3-t2;
tau=tau3-tau1;
P1=cross(rho2_hat,rho3_hat);
P2=cross(rho1_hat,rho3_hat);
P3=cross(rho1_hat,rho2_hat);
D0=dot(rho1_hat,P1);
D11=dot(R1,P1);
D21=dot(R2,P1);
D31=dot(R3,P1);
D12=dot(R1,P2);
D22=dot(R2,P2);
D32=dot(R3,P2);
D13=dot(R1,P3);
D23=dot(R2,P3);
D33=dot(R3,P3);
A=1/D0*(-D12*tau3/tau+D22+D32*tau1/tau);
B=1/6/D0*(D12*(tau3^2-tau^2)*tau3/tau+D32*(tau^2-tau1^2)*tau1/tau);
E=dot(R2,rho2_hat);
R2_pow2=(norm(R2))^2;
a=-(A^2+2*A*E+R2_pow2);
b=-2*muo*B*(A+E);
c=-muo^2*B^2;
syms x;
f(x) = x^8+a*x^6+b*x^3+c;
f_dash(x)=diff(f(x));
sol=f0;
ratio=double(((f(f0))/(f_dash(f0))));
while double(abs(ratio))>= 1e-8;
    sol=double(sol-ratio);
    ratio=double((f(sol))/(f_dash(sol)));
end
r2=sol;
if r2 > 0 && r2 > R
rho1_mag=1/D0*((6*(D31*tau1/tau3+D21*tau/tau3)*r2^3+muo*D31*(tau^2-tau1^2)*tau1/tau3)/(6*r2^3+muo*(tau^2-tau3^2))-D11);
rho2_mag=A+muo*B/r2^3;
rho3_mag=1/D0*((6*(D13*tau3/tau1-D23*tau/tau1)*r2^3+muo*D13*(tau^2-tau3^2)*tau3/tau1)/(6*r2^3+muo*(tau^2-tau1^2))-D33);
r1_vector=R1+rho1_mag*rho1_hat;
r2_vector=R2+rho2_mag*rho2_hat;
r3_vector=R3+rho3_mag*rho3_hat;
f_1=1-1/2*muo/r2^3*tau1^2;
g1=tau1-1/6*muo/r2^3*tau1^3;
f_3=1-1/2*muo/r2^3*tau3^2;
g3=tau3-1/6*muo/r2^3*tau3^3;
v2_vector=1/(f_1*g3-f_3*g1)*(-f_3*r1_vector+f_1*r3_vector);

if improvement == 0
    % orbital element procdure without improvement
    vr2=dot(r2_vector,v2_vector)/norm(r1_vector);
    h=cross(r2_vector,v2_vector);
    mag_h=norm(h);
    i=acosd(h(3)/mag_h);
    N=cross([0,0,1],h);
    mag_N=norm(N);
    if N(2) >= 0
        omega=acosd(N(1)/mag_N);
    elseif N(2) < 0
        omega=360-acosd(N(1)/mag_N);
    end
    e_vector=(((norm(v2_vector))^2-muo/norm(r2_vector))*r2_vector-norm(r2_vector)*vr2*v2_vector)/muo;
    mag_e=norm(e_vector);
    if e_vector(3) >= 0
        w=acosd(dot(N,e_vector)/mag_N/mag_e);
    elseif e_vector(3) < 0
        w=360-acosd(dot(N,e_vector)/mag_N/mag_e);
    end
    if vr2 >= 0
        theta2=acosd(dot(e_vector,r2_vector)/norm(r2_vector)/mag_e);
    elseif vr2 < 0
        theta2=360-acosd(dot(e_vector,r2_vector)/norm(r2_vector)/mag_e);
    end
    rp=mag_h^2/muo/(1+mag_e*cosd(0));
    zp=rp-R;
    epslon=-1/2*muo^2/mag_h^2*(1-mag_e^2);
    
elseif improvement == 1
    rho1_mag_imp_old=0; rho1_mag_imp=1;
    rho2_mag_imp_old=0; rho2_mag_imp=1;
    rho3_mag_imp_old=0; rho3_mag_imp=1;
    while abs(rho1_mag_imp_old-rho1_mag_imp)>1e-8 || abs(rho2_mag_imp_old-rho2_mag_imp)>1e-8 || abs(rho3_mag_imp_old-rho3_mag_imp)>1e-8
        r2_mag_imp=norm(r2_vector);
        v2_mag_imp=norm(v2_vector);
        alpha_imp=2/r2_mag_imp-v2_mag_imp^2/muo;
        vr2_imp=dot(r2_vector,v2_vector)/r2_mag_imp;
        syms n m cay1 cay3;
        z01=0.1; z03=0.1;
        if z01>0 % ellipse
            S1(n)=(sqrt(n)-sin(sqrt(n)))/(sqrt(n))^3;
            C1(m)=(1-cos(sqrt(m)))/m;
        elseif z01<0 % hyperbola
            S1(n)=(sinh(sqrt(-n))-sqrt(-n))/(sqrt(-n))^3;
            C1(m)=(cosh(sqrt(-m))-1)/-m;
        elseif z01==0 % parabola
            S1(n)=sym(1/6);
            C1(m)=sym(1/2);
        end
        if z03>0 % ellipse
            S3(n)=(sqrt(n)-sin(sqrt(n)))/(sqrt(n))^3;
            C3(m)=(1-cos(sqrt(m)))/m;
        elseif z03<0 % hyperbola
            S3(n)=(sinh(sqrt(-n))-sqrt(-n))/(sqrt(-n))^3;
            C3(m)=(cosh(sqrt(-m))-1)/-m;
        elseif z03==0 % parabola
            S3(n)=sym(1/6);
            C3(m)=sym(1/2);
        end
        f1(cay1)=r2_mag_imp*vr2_imp/sqrt(muo)*cay1^2*C1(alpha_imp*cay1^2)+(1-alpha_imp*r2_mag_imp)*cay1^3*S1(alpha_imp*cay1^2)+r2_mag_imp*cay1-sqrt(muo)*tau1;
        f1_dash(cay1)=diff(f1(cay1));
        f3(cay3)=r2_mag_imp*vr2_imp/sqrt(muo)*cay3^2*C3(alpha_imp*cay3^2)+(1-alpha_imp*r2_mag_imp)*cay3^3*S3(alpha_imp*cay3^2)+r2_mag_imp*cay3-sqrt(muo)*tau3;
        f3_dash(cay3)=diff(f3(cay3));
        sol1=z01/alpha_imp;
        ratio1=double(((f1(sol1))/(f1_dash(sol1))));
        while double(abs(ratio1))>= 1e-8;
            sol1=double(sol1-ratio1);
            ratio1=double((f1(sol1))/(f1_dash(sol1)));
            z01=alpha_imp*sol1^2;
            if z01>0 % ellipse
                S1(n)=(sqrt(n)-sin(sqrt(n)))/(sqrt(n))^3;
                C1(m)=(1-cos(sqrt(m)))/m;
            elseif z01<0 % hyperbola
                S1(n)=(sinh(sqrt(-n))-sqrt(-n))/(sqrt(-n))^3;
                C1(m)=(cosh(sqrt(-m))-1)/-m;
            elseif z01==0 % parabola
                S1(n)=sym(1/6);
                C1(m)=sym(1/2);
            end
            f1(cay1)=r2_mag_imp*vr2_imp/sqrt(muo)*cay1^2*C1(alpha_imp*cay1^2)+(1-alpha_imp*r2_mag_imp)*cay1^3*S1(alpha_imp*cay1^2)+r2_mag_imp*cay1-sqrt(muo)*tau1;
            f1_dash(cay1)=diff(f1(cay1));
        end
        Cay1=sol1;
        sol3=z03/alpha_imp;
        ratio3=double(((f3(sol3))/(f3_dash(sol3))));
        while double(abs(ratio3))>= 1e-8;
            sol3=double(sol3-ratio3);
            ratio3=double((f3(sol3))/(f3_dash(sol3)));
            z03=alpha_imp*sol1^2;
            if z03>0 % ellipse
                S3(n)=(sqrt(n)-sin(sqrt(n)))/(sqrt(n))^3;
                C3(m)=(1-cos(sqrt(m)))/m;
            elseif z03<0 % hyperbola
                S3(n)=(sinh(sqrt(-n))-sqrt(-n))/(sqrt(-n))^3;
                C3(m)=(cosh(sqrt(-m))-1)/-m;
            elseif z03==0 % parabola
                S3(n)=sym(1/6);
                C3(m)=sym(1/2);
            end
            f3(cay3)=r2_mag_imp*vr2_imp/sqrt(muo)*cay3^2*C3(alpha_imp*cay3^2)+(1-alpha_imp*r2_mag_imp)*cay3^3*S3(alpha_imp*cay3^2)+r2_mag_imp*cay3-sqrt(muo)*tau3;
            f3_dash(cay3)=diff(f3(cay3));
        end
        Cay3=sol3;
        f1_imp=1-Cay1^2/r2_mag_imp*double(C1(alpha_imp*Cay1^2));
        g1_imp=tau1-1/sqrt(muo)*Cay1^3*double(S1(alpha_imp*Cay1^2));
        f3_imp=1-Cay3^2/r2_mag_imp*double(C3(alpha_imp*Cay3^2));
        g3_imp=tau3-1/sqrt(muo)*Cay3^3*double(S1(alpha_imp*Cay3^2));
        c1=g3_imp/(f1_imp*g3_imp-f3_imp*g1_imp);
        c3=-g1_imp/(f1_imp*g3_imp-f3_imp*g1_imp);
        rho1_mag_imp_old=rho1_mag_imp;
        rho2_mag_imp_old=rho2_mag_imp;
        rho3_mag_imp_old=rho3_mag_imp;
        rho1_mag_imp=1/D0*(-D11+1/c1*D21-c3/c1*D31);
        rho2_mag_imp=1/D0*(-c1*D12+D22-c3*D32);
        rho3_mag_imp=1/D0*(-c1/c3*D13+1/c3*D23-D33);
        r1_vector_imp=R1+rho1_mag_imp*rho1_hat;
        r2_vector_imp=R2+rho2_mag_imp*rho2_hat;
        r3_vector_imp=R3+rho3_mag_imp*rho3_hat;
        v2_vector_imp=1/(f1_imp*g3_imp-f3_imp*g1_imp)*(-f3_imp*r1_vector_imp+f1_imp*r3_vector_imp);
    end
    % orbital element procdure with improvement
    h=cross(r2_vector_imp,v2_vector_imp);
    mag_h=norm(h);
    i=acosd(h(3)/mag_h);
    N=cross([0,0,1],h);
    mag_N=norm(N);
    if N(2) >= 0
        omega=acosd(N(1)/mag_N);
    elseif N(2) < 0
        omega=360-acosd(N(1)/mag_N);
    end
    e_vector=(((norm(v2_vector_imp))^2-muo/norm(r2_vector_imp))*r2_vector_imp-norm(r2_vector_imp)*vr2_imp*v2_vector_imp)/muo;
    mag_e=norm(e_vector);
    if e_vector(3) >= 0
        w=acosd(dot(N,e_vector)/mag_N/mag_e);
    elseif e_vector(3) < 0
        w=360-acosd(dot(N,e_vector)/mag_N/mag_e);
    end
    if vr2_imp >= 0
        theta2=acosd(dot(e_vector,r2_vector_imp)/norm(r2_vector_imp)/mag_e);
    elseif vr2_imp < 0
        theta2=360-acosd(dot(e_vector,r2_vector_imp)/norm(r2_vector_imp)/mag_e);
    end
    rp=mag_h^2/muo/(1+mag_e*cosd(0));
    zp=rp-R;
    epslon=-1/2*muo^2/mag_h^2*(1-mag_e^2);
end
else 
    display('r2 in not acceptable change you initial condition');
end
end