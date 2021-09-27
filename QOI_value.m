function Q_val = QOI_value(lQ)
global lP
global P flag_disp
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2

age_max = P.age_max;
a = P.a;
t = P.t;
da = P.da;
 
FileName = ['Results/SA/',lQ,'_',lP,'_',num2str(P.(lP),5),'.mat'];

if exist(FileName,'file') % Q_val already calculated before
    if flag_disp; disp('load Q_val results...'); end
    load(FileName,'Q_val')
    return
end

switch lQ
    case 'R0'    
        Q_val  = R0_cal();
    case 'RHM'
        [~,Q_val,~]  = R0_cal();
    case 'RMH'
        [~,~,Q_val]  = R0_cal();
    case 'EIR_EE'
        %% EIR
        [SH, EH, DH, AH, SM, EM, IM, ~, ~, ~] = age_structured_Malaria;
        PH_final = SH(:,end)+EH(:,end)+DH(:,end)+AH(:,end); % total human at age a, t = n
        NH = trapz(PH_final)*da;
        NM = SM(end)+EM(end)+IM(end);
        [bh,~] = biting_rate(NH,NM);
        Q_val = bh.*IM(end)./NM*365; % annual EIR
    case 'EE_numerical'
        % use numerical simulation
        [SH, EH, DH, AH, SM, EM, IM, ~, ~, ~] = age_structured_Malaria;
        PH = SH+EH+DH+AH; % total human at age a, t = n
        F_number = @(x) human_model_der_fun2(x);
        x_EE = [SH(:,end);EH(:,end);DH(:,end);AH(:,end)];
        if norm(F_number(x_EE)) > 10^-5
            disp('EE not achieved')
            keyboard
        end
        F_prop = @(x) human_model_der_fun(x);
        x_EE_prop = [SH(:,end)./PH(:,end);EH(:,end)./PH(:,end);DH(:,end)./PH(:,end);AH(:,end)./PH(:,end)];
        norm(F_prop(x_EE_prop))
        keyboard
        %% plot
        figure_setups; hold on;
        plot(a,SH(:,end),'-','Color',colour_mat1);
        plot(a,EH(:,end),'-','Color',colour_mat3);
        plot(a,DH(:,end),'-','Color',colour_mat2);
        plot(a,AH(:,end),'-','Color',colour_mat7);
        plot(a,PH(:,end),'-k');
        legend('SH (numercial)','EH (numercial)','AH (numercial)', 'DH (numercial)','PH (numercial)');
        title('Final Age Dist.');
        xlabel('age');
        axis_years(gca,age_max); % change to x-axis to years if needed
        axis([0 age_max 0 max(SH(:,end)+EH(:,end)+DH(:,end)+AH(:,end))]);
        grid on; grid minor
        figure_setups; hold on;
        plot(a,SH(:,end)./PH(:,end),'-','Color',colour_mat1);
        plot(a,EH(:,end)./PH(:,end),'-','Color',colour_mat3);
        plot(a,DH(:,end)./PH(:,end),'-','Color',colour_mat2);
        plot(a,AH(:,end)./PH(:,end),'-','Color',colour_mat7);
        plot(a,PH(:,end)./PH(:,end),'-k');
        legend('SH prop (numercial)','EH prop (numercial)','AH prop (numercial)', 'DH prop (numercial)','PH (numercial)');
    case 'EE_fsolve'
        F_number = @(x) human_model_der_fun(x);
%         F = @(x) human_model_der_fun2(x);
        x0 = [P.PH_stable(1); 0.1*ones(length(P.a)-1,1).*P.PH_stable(2:end); 0; 0.07*ones(length(P.a)-1,1).*P.PH_stable(2:end);...
            0; 0.2*ones(length(P.a)-1,1).*P.PH_stable(2:end); 0; 0.7*ones(length(P.a)-1,1).*P.PH_stable(2:end)]; 
        error_tolerance = 1e-20;
        options = optimoptions('fsolve','Display','none','OptimalityTolerance', error_tolerance);
        [xsol,error,~,~,~] = fsolve(F_number,x0,options);
        keyboard
        x_EE = reshape(xsol,[P.na,4]);
        %% plot
        figure_setups; hold on;
%         plot(a,x_EE(:,1).*P.PH_stable,'-','Color',colour_mat1); 
%         plot(a,x_EE(:,2).*P.PH_stable,'-','Color',colour_mat3);
%         plot(a,x_EE(:,3).*P.PH_stable,'-','Color',colour_mat2);
%         plot(a,x_EE(:,4).*P.PH_stable,'-','Color',colour_mat7);
%         plot(a,sum(x_EE,2).*P.PH_stable,'-k');
%         axis([0 age_max 0 max(sum(x_EE,2).*P.PH_stable)]);
        plot(a,x_EE(:,1),'-','Color',colour_mat1); 
        plot(a,x_EE(:,2),'-','Color',colour_mat3);
        plot(a,x_EE(:,3),'-','Color',colour_mat2);
        plot(a,x_EE(:,4),'-','Color',colour_mat7);
        plot(a,sum(x_EE,2),'-k');
        keyboard
        axis([0 age_max 0 max(sum(x_EE,2))]);
        legend('SH (solver)','EH (solver)','DH (solver)', 'AH (solver)', 'PH (solver)');
        title('Final Age Dist.');
        xlabel('age');
        axis_years(gca,P.age_max); % change to x-axis to years if needed
        grid on
        figure_setups; hold on;
%         plot(a,x_EE(:,1),'-','Color',colour_mat1); 
%         plot(a,x_EE(:,2),'-','Color',colour_mat3);
%         plot(a,x_EE(:,3),'-','Color',colour_mat2);
%         plot(a,x_EE(:,4),'-','Color',colour_mat7);
%         plot(a,sum(x_EE,2),'-k');

        plot(a,x_EE(:,1)./P.PH_stable,'-','Color',colour_mat1); 
        plot(a,x_EE(:,2)./P.PH_stable,'-','Color',colour_mat3);
        plot(a,x_EE(:,3)./P.PH_stable,'-','Color',colour_mat2);
        plot(a,x_EE(:,4)./P.PH_stable,'-','Color',colour_mat7);
        plot(a,sum(x_EE,2)./P.PH_stable,'-k');

        legend('SH prop (solver)','EH prop (solver)','DH prop (solver)', 'AH prop (solver)', 'PH prop (solver)');
        keyboard
    case 'stability'
        
    otherwise
        keyboard
end

% save(FileName,'Q_val')

end
