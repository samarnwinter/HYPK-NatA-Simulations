% Stochastic model for HYPK–NatA-mediated cotranslational N-terminal acetylation
% -------------------------------------------------------------------------
% This MATLAB script is a Gillespie algorithm-based simulation model of 
% cotranslational N-terminal acetylation (NTA). The system is modeled 
% using discrete reaction events governed by rate constants for NatA–ribosome 
% association (k12, k43), dissociation (k21, k34), acetylation (k23),
% translation elongation rate (ω) and coding sequence length (L).
%------------------------------------------------------------------------
% Usage:
%   Run directly in MATLAB:
%       HYPK_NatA_Ac_sim_model.m
%
% Author: Inayat Ulah Irshad
% Date: October 2025


clear all;clc; close all;
tic
n=300; 
nat=[10 20 40];   
a=1;R_total=70e-9;N_total=nat(1)*(10^-9);
k1=4*1e8;k2=0.1;k4=k2;k5=k1;          %======= hypk-
%k1=0.2*1e8;k2=1;k4=k2;k5=k1;         %======= hypk+
kd=k2/k1;
b=(R_total-N_total+kd); 
P1=[];P2=[];P3=[];P4=[];prnc=[];pnatarnc=[];p3_arnc=[];p4arnc=[];total_acetyl=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%rate calculations  

c=-kd*N_total;

d=b^2-4*a*c;
if d>0
    NatA_s1=(-b+sqrt(d))/2*a;
    NatA_s2=(-b-sqrt(d))/2*a;
end


NC_length = [68, 88, 95, 105, 115, 135, 155, 185, 205];
Avg_k = [0, 0.001373, 0.002281, 0.845, 0.2825, 0.17, 0.285, 0.075, 0.046];
a1 = 0.79; a2 = 0.2; a3 = 0.2950;
b1 = 105.0; b2 = 116.2; b3 = 154.5404;
c1 = 6.414; c2 = 10.6473; c3 = 28.1588;
x = 1:1:353;
y = a1 * exp(-((x - b1) / c1).^2) + a2 * exp(-((x - b2) / c2).^2) + a3 * exp(-((x - b3) / c3).^2);
y=y/max(y);
%plot(x,y); hold on; plot(NC_length,Avg_k);

kk1=(k1)*NatA_s1;kk2=k2;kk3=1.0;kk4=kk2;kk5=kk1;



switchx="on";                             % to start estimating time of simulation

no_of_iteration = 5000;

                                          % length of mrna track of a particular gen


W=zeros(1,n);
W(1:end)=4;
alpha_rate=0.4*ones(1,100);
m=length(alpha_rate);                     % for uniform translational rate
p=length(W);                              % length of mrna track
w= W;                                     
L=10;                                     % length of ribosome

flux=zeros(1,p);                          % Count total number of ribosomes on each codon for all values of alpha, therefore initailizing at each codon as 0
flux2=zeros(1,p);                         % Count number of ribosome on each codon for a particular value of alpha, therefore initailizing at each codon as 0
flux3=zeros(1,m);                         % Count number of ribosomes at 2nd codon for each value of alpha,therefore initailizing counts at 2nd codon as 0 for each alpha value
flux5=zeros(1,m);                         % Count number of ribosomes at pth codon for each value of alpha,therefore initailizing counts at 2nd codon as 0 for each alpha value
reads1=zeros(1,p);                        % Initializing reads on each of the codon
reads2=zeros(1,p);                        % Initializing reads on each of the codon
reads3=zeros(1,p);                        % Initializing reads on each of the codon
reads4=zeros(1,p);                        % Initializing reads on each of the codon

avg_density=zeros(1,m); 
total_flux_termination=zeros(1,m);           
total_flux_initiation=zeros(1,m);

t= zeros(1,m);
state_ac(1:p)={[0 0 0 0]};

for i=1:n
    kk3=max(Avg_k);
    rate_ac{i}=[kk1,kk2,kk3*y(i),kk4,kk5];
    rate_k3(i)=kk3*y(i);
end 
rate_ac(1);
rr=0;
rr1=0;
q=1;
xm=1;
time_prnc=[];time_pnatarnc=[];time_p3_arnc=[];time_p4arnc=[];time_total_acetyl=[];time_time=[];

check="nexecuted";
figure(1); hold on;

for k=1:100                                 % loop over alpha(initiation rates)
    state_ac(1:p+L)={[0 0 0 0]};
    flux = zeros(1,p);                     % initializing total counts on each codon as 0 for each value of initiation rate
    state = zeros(1,p+L);                  % represents the mrna stand
    time=0;
    s=100*rand();                          % generating a random number between 0 to 100
    counter=0;
    alpha_rate=0.4*ones(1,100);
    k
    cm=0;
    clear time_ac;
    clear act_st;
    for i=1:no_of_iteration               % moving the loop on same mrna track this many times
        if counter==1
            alpha_rate=zeros(1,100);
        end    
        if sum(state_ac{p})==1
            break
        end    
   
        
        
        s=100*rand();                     % generating a random number between 0 to 100
        beta= 9.0;                        % termination rate
        i_ribo= 0;                        % ribosome indexing ( position of ribosome on mrna track)
        R = alpha_rate(k);                % initiation rate
        t_occup=0;     
        
        for j = 2:p-1
            if sum(state_ac{j})>0 && sum(cell2mat(state_ac(j+1:j+L)))==0
                t_occup = t_occup + 1;
                i_ribo(t_occup) = j;

                if state_ac{j}(1)==1
                    R = R + rate_ac{j}(1) + w(j);
                elseif state_ac{j}(2)==1
                    R = R + rate_ac{j}(2) + rate_ac{j}(3) + w(j);
                elseif state_ac{j}(3)==1
                    R = R + rate_ac{j}(4) + w(j);
                elseif state_ac{j}(4)==1
                    R = R + rate_ac{j}(5) + w(j);
                end
            end
        end

        if sum(state_ac{p}(1:4))==1                   % for last site p if it is occupied then
            R= R + beta ;                             % adding the termination rate
            t_occup=t_occup+1;                        % again add the ribosomes
            i_ribo(t_occup)=p;                        % and do the indexing
            
        end
        
        r1= R*rand();                           % random no. generation b/w 0 & R 
        r2=rand();                              % random no. b/w 0 and 1
        tau=(1/R)*log(1/r2);                    % exponential distribution time
        time=time+tau;     
        
        initiation=0;
        rate_total=alpha_rate(k);
        if r1>0 && r1<=alpha_rate(k)
            
            counter=counter+1;
            for m=2:11
                initiation=initiation+sum(state_ac{m}(1:4));
            end
            if initiation==0
                state_ac{2}(1)=1;
                flux5(k)=flux5(k)+1;
            end 
            
            
        else
            for j=1:t_occup
                rate=rate_total;
                if i_ribo(j)==p 
                    
                    if r1>rate && r1<=rate+beta
                        state_ac{i_ribo(j)}(1:4)=0;
                        flux3(k)=flux3(k)+1;
                       
                    end
                    rate=rate+beta;



                elseif state_ac{i_ribo(j)}(1)==1
                    if r1>rate && r1<= rate+w(i_ribo(j))

                        if (state_ac{i_ribo(j)}(1)==1 && sum(state_ac{i_ribo(j)+L}))==0 
                            state_ac{i_ribo(j)+1}=state_ac{i_ribo(j)};
                            state_ac{i_ribo(j)}(1:4)=0;
                            %break;
                            
                        end
                    
                        
                    elseif r1>rate+w(i_ribo(j)) && r1<= rate+w(i_ribo(j))+rate_ac{i_ribo(j)}(1)
                        state_ac{i_ribo(j)}(1:4)=0;
                        state_ac{i_ribo(j)}(2)=1;
                        break;
                    end
                    
                    rate=rate+w(i_ribo(j))+rate_ac{i_ribo(j)}(1);

                elseif state_ac{i_ribo(j)}(2)==1

                    if r1>rate && r1<= rate+w(i_ribo(j))

                        if (state_ac{i_ribo(j)}(2)==1 && sum(state_ac{i_ribo(j)+L}))==0 
                            state_ac{i_ribo(j)+1}(1:4)=state_ac{i_ribo(j)}(1:4);
                            state_ac{i_ribo(j)}(1:4)=0;
                            %break;
                            
                        end
                    elseif r1>rate+w(i_ribo(j)) && r1<= rate+w(i_ribo(j))+rate_ac{i_ribo(j)}(2)
                        state_ac{i_ribo(j)}(1:4)=0;
                        state_ac{i_ribo(j)}(1)=1;
                        break;
                   
                    elseif r1>rate+w(i_ribo(j))+rate_ac{i_ribo(j)}(2) && r1<= rate+w(i_ribo(j))+rate_ac{i_ribo(j)}(2)+rate_ac{i_ribo(j)}(3)
                        state_ac{i_ribo(j)}(1:4)=0;
                        state_ac{i_ribo(j)}(3)=1;  
                        break;

                    end 
                    rate=rate+w(i_ribo(j))+rate_ac{i_ribo(j)}(2)+rate_ac{i_ribo(j)}(3);
                
                elseif state_ac{i_ribo(j)}(3)==1   
                    
                    if r1>rate && r1<= rate+w(i_ribo(j))

                        if (state_ac{i_ribo(j)}(3)==1 && sum(state_ac{i_ribo(j)+L}))==0 
                            state_ac{i_ribo(j)+1}(1:4)=state_ac{i_ribo(j)}(1:4);
                            state_ac{i_ribo(j)}(1:4)=0;
                            break
                            
                        end
                    elseif r1>rate+w(i_ribo(j)) && r1<= rate+w(i_ribo(j))+rate_ac{i_ribo(j)}(4)
                        state_ac{i_ribo(j)}(1:4)=0;
                        state_ac{i_ribo(j)}(4)=1;
                        break
            
                    end 
                    rate=rate+w(i_ribo(j))+rate_ac{i_ribo(j)}(4);
               
                elseif state_ac{i_ribo(j)}(4)==1   
                    
                    if r1>rate && r1<= rate+w(i_ribo(j))

                        if (state_ac{i_ribo(j)}(4)==1 && sum(state_ac{i_ribo(j)+L}))==0 
                            state_ac{i_ribo(j)+1}(1:4)=state_ac{i_ribo(j)}(1:4);
                            state_ac{i_ribo(j)}(1:4)=0;
                            break;
                            
                        end
                    elseif r1>rate+w(i_ribo(j)) && r1<= rate+w(i_ribo(j))+rate_ac{i_ribo(j)}(5)
                        state_ac{i_ribo(j)}(1:4)=0;
                        state_ac{i_ribo(j)}(3)=1;
                        break
            
                    end 
                    rate=rate+w(i_ribo(j))+rate_ac{i_ribo(j)}(5);
                      
                
                    
                end
                rate_total=rate;   
            end 
        end    
          
        total_flux_termination(k)=flux3(k)/time;
        total_flux_initiation(k)=flux5(k)/time;


       %%%%record time trajectories

        
        cm=cm+1; 
        index = find(cellfun(@(x) sum(x) == 1, state_ac));
        %index
        time_ac(cm)=time;
        if state_ac{index}(1)==1
            act_st(cm)=0;
        elseif state_ac{index}(2)==1
            act_st(cm)=1;
        elseif state_ac{index}(3)+state_ac{index}(4)==1
            act_st(cm)=2;
        end  
    end
    %plot(time_ac,act_st, 'LineWidth', 1.5);   hold on 
    tim_d{k}=time_ac;
    act_d{k}=act_st;
       
end    

hold off;
j1=alpha_rate.*(1-alpha_rate./5);
j2=(1+alpha_rate.*9/5);
j=j1./j2;
reads=reads1+reads2+reads3+reads4;
prnc=reads1(2:end)./reads(2:end);
pnatarnc=reads2(2:end)./reads(2:end);
p3_arnc=reads3(2:end)./reads(2:end);
p4arnc=reads4(2:end)./reads(2:end);
total_acetyl=p3_arnc+p4arnc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

kmm=0;kmm1=0;
for i=1:100
    
    if any(act_d{i} > 1);
        kmm=kmm+1;
        st(kmm)=i;
    else
        kmm1=kmm1+1;
        st1(kmm1)=i;
    end   
end



%%
figure(1); 
hold on;
legendEntries = cell(1, 10); % Preallocate legend entries

% Define 10 distinct colors
colors = [0 0 1; 1 0 0; 0 0.6 0; 0.9 0.6 0; 0.5 0 0.5; 
          0.2 0.8 0.8; 0.8 0 0.8; 0.6 0.6 0; 0.3 0.3 1; 1 0.3 0.35];

for i = 1:10
    if rand() >= 0.5
        pt = st(randi(length(st)));
    else    
        pt = st1(randi(length(st1)));
    end

    plot(tim_d{pt}, act_d{pt}, '--', 'LineWidth', 1.0, 'Color', colors(i, :)); 
    legendEntries{i} = sprintf('T%d', i); % Store legend label
end    

hold off;
xlabel('Time (s)', 'FontSize', 16, 'FontName', 'Arial');
yticks([0 1 2]); % Assign positions (adjust as needed)
yticklabels({'RNC', 'RNC\cdot NatA', 'Ac-RNC\cdotNatA'}); 
%ytickangle(90);

set(gca, 'TickLabelInterpreter', 'tex');
set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.015 0.015], 'FontName', 'Arial');
leg = legend(legendEntries, 'Location', 'best', 'FontSize', 14, 'Interpreter', 'tex');
set(leg, 'Position', [0.8, 0.5, 0.2, 0.2]); % Adjust values as needed
set(leg, 'Box', 'off', 'FontName', 'Arial');
hold off;


