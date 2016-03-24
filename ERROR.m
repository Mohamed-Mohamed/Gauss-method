function [ e,Max_e,std_e, mean_e, RMS_e ] = ERROR ( a1,a2 )
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg

e=a1-a2;
Max_e=max(abs(e));
std_e=std(e);
mean_e=mean(e);
RMS_e=rms(e);
end

