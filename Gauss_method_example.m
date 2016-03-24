%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg



% Gauss method example
close all; clear all; clc;
%% inputs
rho1_hat=[0.71643, 0.68074, -0.1527];               % direction cosine vectors of observer 1 (1x3)
rho2_hat=[0.56897, 0.79531, -0.20917];               % direction cosine vectors of observer 2 (1x3)
rho3_hat=[0.41841, 0.87007, -0.26059];               % direction cosine vectors of observer 3 (1x3)
R1=[3489.8, 3430.2, 4078.5];                         % position of observer 1 (1x3)
R2=[3460.1, 3460.1, 4078.5];                         % position of observer 2 (1x3)
R3=[3429.9, 3490.1, 4078.5];                         % position of observer 3 (1x3)
t1=0;                          % time of observer 1 in sec
t2=118.1;                          % time of observer 2 in sec
t3=237.58;                          % time of observer 3 in sec
f0=9e3;                          % initial value of r2
noprovement=0;       % if 0 no improve in orbital elements, if 1 improve in orbital elements
improvement=1;       % if 0 no improve in orbital elements, if 1 improve in orbital elements
muo=398600; %  Gravitational Parameter
Re=6378;    % raduis of plant
%% Gauss without improvement Solution
[ h, mag_h, i, omega, e_vector, mag_e, w, theta, rp, zp, epslon ] = Gauss ( rho1_hat, rho2_hat, rho3_hat, R1, R2, R3, t1, t2, t3, muo, f0, Re, noprovement );
%% Gauss with improvement Solution
[ h_imp, mag_h_imp, i_imp, omega_imp, e_vector_imp, mag_e_imp, w_imp, theta_imp, rp_imp, zp_imp, epslon_imp ] = Gauss ( rho1_hat, rho2_hat, rho3_hat, R1, R2, R3, t1, t2, t3, muo, f0, Re, improvement );
%% orbital elements
Gauss_Without_Improvment=[mag_h; i; omega; mag_e; w; theta];
Gauss_With_Improvment=[mag_h_imp; i_imp; omega_imp; mag_e_imp; w_imp; theta_imp];
f = figure(8);
set(gcf,'color','w');
set(f,'Position',[100 100 450 250])
% create the data
d = [Gauss_Without_Improvment, Gauss_With_Improvment];
% Create the column and row names in cell arrays 
cnames = {'Gauss Without Improvment','Gauss With Improvment'};
rnames = {'h','i','Omega','e','w','theta'};
% Create the uitable
t = uitable(f,'Data',d,...
            'ColumnName',cnames,... 
            'RowName',rnames,'ColumnWidth',{170},'Position',[20 100 408 130]);
%% RK4 parameter of orbital element without improvement Solution
order=6;
X0=[mag_h;mag_e;i;omega;w;theta];
B=[0;0;0;0;0;0;];
sol(1:6,1)=X0;
dt=10000;
t_initial=0;
t_final=3e6+4e5*0;
%% RK4 parameter of orbital element with improvement Solution
X0_imp=[mag_h_imp;mag_e_imp;i_imp;omega_imp;w_imp;theta_imp];
B_imp=[0;0;0;0;0;0;];
sol_imp(1:6,1)=X0_imp;
%% solution of orbital element without improvement Solution
for n=1:length(t_initial:dt:t_final)
    A=[0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      muo^2/sol(1,n)^4*(1+sol(2,n)*cosd(sol(6,n)))^2,0,0,0,0,0];    
    [ XX, t ] = RK4( A,B,sol(1:6,n),dt,n*dt,(n+1)*dt,order );
    sol(1:6,n+1)=XX(1:6,2);
end
for m=1:length(sol(1,:))
    [ r, v ] = OrbitalElements2rvGeo( sol(1,m), muo, sol(2,m), sol(3,m), sol(4,m), sol(5,m), sol(6,m) );
    r_XYZ(1:3,m)=r;
    v_XYZ(1:3,m)=v;
end
%% solution of orbital element with improvement Solution
for n=1:length(t_initial:dt:t_final)
    A_imp=[0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      muo^2/sol_imp(1,n)^4*(1+sol_imp(2,n)*cosd(sol_imp(6,n)))^2,0,0,0,0,0];    
    [ XX_imp, t_imp ] = RK4( A_imp,B_imp,sol_imp(1:6,n),dt,n*dt,(n+1)*dt,order );
    sol_imp(1:6,n+1)=XX_imp(1:6,2);
end
for m=1:length(sol_imp(1,:))
    [ r_imp, v_imp ] = OrbitalElements2rvGeo( sol_imp(1,m), muo, sol_imp(2,m), sol_imp(3,m), sol_imp(4,m), sol_imp(5,m), sol_imp(6,m) );
    r_XYZ_imp(1:3,m)=r_imp;
    v_XYZ_imp(1:3,m)=v_imp;
end
%% RK4 parameter of orbital element with oblateness without improvement Solution
X01=[mag_h;mag_e;i;omega;w;theta];
B1=[0;0;0;0;0;0;];
sol1(1:6,1)=X0;
J2=1.08263e-3;
%% RK4 parameter of orbital element with oblateness with improvement Solution
X01_imp=[mag_h_imp;mag_e_imp;i_imp;omega_imp;w_imp;theta_imp];
B1_imp=[0;0;0;0;0;0;];
sol1_imp(1:6,1)=X0_imp;
%% solution of orbital element with oblateness without improvement Solution
for n=1:length(t_initial:dt:t_final)
    A1=[0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      -(3/2/sol1(1,n)^8*muo^4*J2*Re^2*(1-sol1(2,n)^2)^(3/2))*cosd(sol1(3,n)),0,0,0,0,0;...
      -(3/2/sol1(1,n)^8*muo^4*J2*Re^2*(1-sol1(2,n)^2)^(3/2))*(5/2*(sind(sol1(3,n)))^2-2),0,0,0,0,0;...
      muo^2/sol1(1,n)^4*(1+sol1(2,n)*cosd(sol1(6,n)))^2,0,0,0,0,0];    
    [ XX1, t1 ] = RK4( A1,B1,sol1(1:6,n),dt,n*dt,(n+1)*dt,order );
    sol1(1:6,n+1)=XX1(1:6,2);
end
for m=1:length(sol1(1,:))
    [ r1, v1 ] = OrbitalElements2rvGeo( sol1(1,m), muo, sol1(2,m), sol1(3,m), sol1(4,m), sol1(5,m), sol1(6,m) );
    r_XYZ1(1:3,m)=r1;
    v_XYZ1(1:3,m)=v1;
end
%% solution of orbital element with oblateness with improvement Solution
for n=1:length(t_initial:dt:t_final)
    A1_imp=[0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      0,0,0,0,0,0;...
      -(3/2/sol1_imp(1,n)^8*muo^4*J2*Re^2*(1-sol1_imp(2,n)^2)^(3/2))*cosd(sol1_imp(3,n)),0,0,0,0,0;...
      -(3/2/sol1_imp(1,n)^8*muo^4*J2*Re^2*(1-sol1_imp(2,n)^2)^(3/2))*(5/2*(sind(sol1_imp(3,n)))^2-2),0,0,0,0,0;...
      muo^2/sol1_imp(1,n)^4*(1+sol1_imp(2,n)*cosd(sol1_imp(6,n)))^2,0,0,0,0,0];    
    [ XX1_imp, t1_imp ] = RK4( A1_imp,B1_imp,sol1_imp(1:6,n),dt,n*dt,(n+1)*dt,order );
    sol1_imp(1:6,n+1)=XX1_imp(1:6,2);
end
for m=1:length(sol1_imp(1,:))
    [ r1_imp, v1_imp ] = OrbitalElements2rvGeo( sol1_imp(1,m), muo, sol1_imp(2,m), sol1_imp(3,m), sol1_imp(4,m), sol1_imp(5,m), sol1_imp(6,m) );
    r_XYZ1_imp(1:3,m)=r1_imp;
    v_XYZ1_imp(1:3,m)=v1_imp;
end
%% plotting
figure(1);
view(3);
set(gcf,'color','w');
subplot(1,2,1)
view(3)
plot3(r_XYZ(1,:),r_XYZ(2,:),r_XYZ(3,:),'linewidth',2);
grid on;
xlabel('X','fontsize',18);
ylabel('Y','fontsize',18);
zlabel('Z','fontsize',18);
title('Solution without oblateness and without improvement','fontsize',18);
subplot(1,2,2)
plot3(r_XYZ1(1,:),r_XYZ1(2,:),r_XYZ1(3,:),'linewidth',2);
grid on;
xlabel('X','fontsize',18);
ylabel('Y','fontsize',18);
zlabel('Z','fontsize',18);
title('Solution with oblateness and without improvement','fontsize',18);

figure(3);
view(3);
set(gcf,'color','w');
subplot(1,2,1)
plot3(r_XYZ_imp(1,:),r_XYZ_imp(2,:),r_XYZ_imp(3,:),'linewidth',2);
grid on;
xlabel('X','fontsize',18);
ylabel('Y','fontsize',18);
zlabel('Z','fontsize',18);
title('Solution without oblateness and with improvement','fontsize',18);
subplot(1,2,2)
plot3(r_XYZ1_imp(1,:),r_XYZ1_imp(2,:),r_XYZ1_imp(3,:),'linewidth',2);
grid on;
xlabel('X','fontsize',18);
ylabel('Y','fontsize',18);
zlabel('Z','fontsize',18);
title('Solution with oblateness and with improvement','fontsize',18);

figure(5);
set(gcf,'color','w');
subplot(1,2,1)
hold all;
view(3);
plot3(r_XYZ(1,:),r_XYZ(2,:),r_XYZ(3,:),'linewidth',2);
plot3(r_XYZ_imp(1,:),r_XYZ_imp(2,:),r_XYZ_imp(3,:),'r','linewidth',2);
grid on;
xlabel('X','fontsize',18);
ylabel('Y','fontsize',18);
zlabel('Z','fontsize',18);
title('Solution without oblateness','fontsize',18);
legend('solution without improvement','solution with improvement');
subplot(1,2,2)
hold all;
view(3);
plot3(r_XYZ1(1,:),r_XYZ1(2,:),r_XYZ1(3,:),'linewidth',2);
plot3(r_XYZ1_imp(1,:),r_XYZ1_imp(2,:),r_XYZ1_imp(3,:),'r','linewidth',2);
grid on;
xlabel('X','fontsize',18);
ylabel('Y','fontsize',18);
zlabel('Z','fontsize',18);
title('Solution with oblateness','fontsize',18);
legend('solution without improvement','solution with improvement');
%% error 
for k=1:length(r_XYZ1(2,:))
    E1(1,k)=norm([r_XYZ(1:3,k)]);
    E2(1,k)=norm([r_XYZ1(1:3,k)]);
end
[ E,Max_e,std_e, mean_e, RMS_e ] = ERROR ( E1(1,:),E2 );
figure(2);
set(gcf,'color','w');
plot(t_initial:dt:t_final+dt,E)
xlim([t_initial, t_final+dt])
xlabel('Time (sec)','fontsize',18);
ylabel('Sol_w_i_t_h_o_u_t _o_b_l-Sol_w_i_t_h _o_b_l','fontsize',18);
title('Error without improvement','fontsize',18);
legend(['Max.(Error) = ' num2str(Max_e) ' Std(Error) = ' num2str(std_e) ' mean(Error) = ' num2str(mean_e) ' RMS(Error) = ' num2str(RMS_e) ]);
grid on;

for k=1:length(r_XYZ1_imp(2,:))
    E1_imp(1,k)=norm([r_XYZ_imp(1:3,k)]);
    E2_imp(1,k)=norm([r_XYZ1_imp(1:3,k)]);
end
[ E_imp,Max_e_imp,std_e_imp, mean_e_imp, RMS_e_imp ] = ERROR ( E1_imp,E2_imp );
figure(4);
set(gcf,'color','w');
plot(t_initial:dt:t_final+dt,E_imp)
xlim([t_initial, t_final+dt])
xlabel('Time (sec)','fontsize',18);
ylabel('Sol_w_i_t_h_o_u_t _o_b_l-Sol_w_i_t_h _o_b_l','fontsize',18);
title('Error with improvement','fontsize',18);
legend(['Max.(Error) = ' num2str(Max_e_imp) ' Std(Error) = ' num2str(std_e_imp) ' mean(Error) = ' num2str(mean_e_imp) ' RMS(Error) = ' num2str(RMS_e_imp) ]);
grid on;

for k=1:length(r_XYZ1_imp(2,:))
    E11(1,k)=norm([r_XYZ(1:3,k)]);
    E22(1,k)=norm([r_XYZ_imp(1:3,k)]);
end
[ E12,Max_e_12,std_e_12, mean_e_12, RMS_e_12 ] = ERROR ( E11,E22 );
figure(6);
set(gcf,'color','w');
plot(t_initial:dt:t_final+dt,E12)
xlim([t_initial, t_final+dt])
xlabel('Time (sec)','fontsize',18);
ylabel('Sol_w_i_t_h_o_u_t _i_m_p-Sol_w_i_t_h _i_m_p','fontsize',18);
title('Error without oblateness','fontsize',18);
legend(['Max.(Error) = ' num2str(Max_e_12) ' Std(Error) = ' num2str(std_e_12) ' mean(Error) = ' num2str(mean_e_12) ' RMS(Error) = ' num2str(RMS_e_12) ]);
grid on;

for k=1:length(r_XYZ1_imp(2,:))
    E11_imp(1,k)=norm([r_XYZ1(1:3,k)]);
    E22_imp(1,k)=norm([r_XYZ1_imp(1:3,k)]);
end
[ E12_imp,Max_e_12_imp,std_e_12_imp, mean_e_12_imp, RMS_e_12_imp ] = ERROR ( E11_imp,E22_imp );
figure(7);
set(gcf,'color','w');
plot(t_initial:dt:t_final+dt,E12_imp)
xlim([t_initial, t_final+dt])
xlabel('Time (sec)','fontsize',18);
ylabel('Sol_w_i_t_h_o_u_t _i_m_p-Sol_w_i_t_h _i_m_p','fontsize',18);
title('Error with oblateness','fontsize',18);
legend(['Max.(Error) = ' num2str(Max_e_12_imp) ' Std(Error) = ' num2str(std_e_12_imp) ' mean(Error) = ' num2str(mean_e_12_imp) ' RMS(Error) = ' num2str(RMS_e_12_imp) ]);
grid on;