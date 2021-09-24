global P
global colour_mat1 colour_mat2 colour_mat3 colour_mat4 colour_mat5 colour_mat6 colour_mat7
global colour_r1 colour_r2
tic
%% setup parameters, stable age dist, etc.
tfinal = 100*365; % final time in days
age_max = 60*365; % max ages in days
P.age_max = age_max;
dt = 20; % time/age step size in days, default = 5;
da = dt;
t = (0:dt:tfinal)';
nt = length(t);
a = (0:da:age_max)';
na = length(a);

P.a = a;
P.na = na;
P.nt = nt;
P.dt = dt;
P.da = da;

% model parameters - rates are in 1/day
Malaria_parameters_baseline;
%% Set range for the bifurcation parameter
param = 0.1125:0.025:0.25; % changing P.betaM the mosquito desired bites
%% options 
error_tolerance = 1e-20;
options = optimoptions('fsolve','Display','none','OptimalityTolerance',...
    error_tolerance);
plot_equilibrium = 1; % can set to zero if working with the DFE
%% Check if DFE is stable numerically
re_max_endemic = zeros(1,length(param));
sol_norm_endemic = zeros(1,length(param));
error_endemic = zeros(1,length(param));
re_max_DFE = zeros(1,length(param));
sol_norm_DFE = zeros(1,length(param));
error_DFE = zeros(1,length(param));
for i = 1:length(param)
    disp(['progress = ',num2str((i-1)/(length(param))*100),'%']);
    P.betaM = param(i);
    P.phi = 1/2;
    P.psi = 1/2;
    P.rho = 1/2;
    [bH,bM] = biting_rate(1,P.gM/P.muM);
    Lambda = 1/(da.*trapz(exp(-P.p_hat*P.a-P.muH_int)));
    Lambda_M = @(x) bM*Lambda*da*trapz(exp(-P.muH_int).*(P.betaD*x(:,3)+P.betaA*x(:,4)));
    Lambda_H = @(x) bH*P.betaM*(P.sigma/(P.sigma+P.muM))*((Lambda_M(x))/(Lambda_M(x) + P.muM));
    F_PSI = @(x) [[-x(1,1)+1; -Lambda_H(x)*x(2:end,1) + P.phi*P.rD*x(2:end,3) + P.rA*x(2:end,4) - diff(x(:,1))/da];...
        [-x(1,2); Lambda_H(x)*x(2:end,1) - P.h*x(2:end,2) - diff(x(:,2))/da];...
        [-x(1,3); P.rho*P.h*x(2:end,2) + P.psi*Lambda_H(x)*x(2:end,4) - P.rD*x(2:end,3) - diff(x(:,3))/da];...
        [-x(1,4); (1-P.rho)*P.h*x(2:end,2) - P.psi*Lambda_H(x)*x(2:end,4) + (1-P.phi)*P.rD*x(2:end,3) - P.rA*x(2:end,4) - diff(x(:,4))/da]];
    % solve for the endemic equilibrium
    P1 = [0.4*ones(length(a),1), 0.2*ones(length(a),1), 0.2*ones(length(a),1), 0.2*ones(length(a),1)]; % initial guess for the DFE
    [eq_age,fval,~,~,jacobian] = fsolve(F_PSI,P1,options);    
    sol_norm_endemic(i) = da*trapz(eq_age(:,1))/365; % convert to years
    if plot_equilibrium == 1
        % plot the solution of the nonlinear solver
        figure_setups;
        plot(a,eq_age(:,1),'-.','Color',colour_mat1); hold on;
        plot(a,eq_age(:,2),'-.','Color',colour_mat3);
        plot(a,eq_age(:,4),'-.','Color',colour_mat2);
        plot(a,eq_age(:,3),'-.','Color',colour_mat7);
        plot(a,sum(eq_age,2),'-.k');
        legend('SH (solver)','EH (solver)','AH (solver)', 'DH (solver)','PH (solver)');
        title('Final Age Dist. Proportions');
        xlabel('age');
        axis_years(gca,age_max); % change to x-axis to years if needed
        grid on
        axis([0 age_max 0 max(sum(eq_age,2))*1.1]);
    end
    error_endemic(i) = max(max(fval));
    temp_eig = sort(real(eig(jacobian)),'descend');
    re_max_endemic(i) = max(temp_eig);
    %if re_max_endemic(i) < 0
    %    disp(['max real part of eigenvalues at DFE = ',num2str(re_max_endemic(i),'%10.6f'), '; endemic equilibrium is stable']);
    %else
    %    disp(['max real part of eigenvalues = ',num2str(re_max_endemic(i),'%10.6f'), '; endemic equilibrium is unstable']);
    %end
    % analyze the DFE
    P2 = [ones(length(a),1), 0*ones(length(a),1), 0*ones(length(a),1), 0*ones(length(a),1)]; % initial guess for the DFE
    [eq_age,fval,~,~,jacobian] = fsolve(F_PSI,P2,options);
    sol_norm_DFE(i) = da*trapz(eq_age(:,1))/365; % convert to years
    error_DFE(i) = max(max(fval));
    temp_eig = sort(real(eig(jacobian)),'descend');
    re_max_DFE(i) = max(temp_eig);
end
%figure_setups;
% scatter(real(eig(jacobian)),imag(eig(jacobian)));
% grid on;
% xline(0,'r','LineWidth',2);
% title('Linearized spectrum at DFE');
% xlim([min(temp_eig) re_max+0.1]);
%% Plot the results
%figure_setups;
figure(9);
hold on;
error_tolerance = 1e-6;
% endemic equilibrium plotting
for j = 1:length(param)
    if re_max_endemic(j) < 0 && error_endemic(j) < error_tolerance
        scatter(param(j),sol_norm_endemic(j),'b','filled');
    elseif re_max_endemic(j) < 0 && error_endemic(j) > error_tolerance
        scatter(param(j),sol_norm_endemic(j),'r','filled');
    elseif re_max_endemic(j) > 0 && error_endemic(j) > error_tolerance
        scatter(param(j),sol_norm_endemic(j),'r');
    else
        scatter(param(j),sol_norm_endemic(j),'b');
    end
end
% DFE plotting
for j = 1:length(param)
    if re_max_DFE(j) < 0 && error_DFE(j) < error_tolerance
        scatter(param(j),sol_norm_endemic(j),'b','filled');
    elseif re_max_DFE(j) < 0 && error_DFE(j) > error_tolerance
        scatter(param(j),sol_norm_DFE(j),'r','filled');
    elseif re_max_DFE(j) > 0 && error_DFE(j) > error_tolerance
        scatter(param(j),sol_norm_DFE(j),'r');
    else
        scatter(param(j),sol_norm_DFE(j),'b');
    end
end
grid on;
xlabel('$\beta_M$');
ylabel('$|| \tilde{S}_H ||$');
set(get(gca,'ylabel'),'rotation',0)
%%
toc