function  [SH, EH, DH, AH, SM, EM, IM, Cm, Cac, Ctot] = age_structured_Malaria_IC(state)
% define the initial condition of the simulation

global P

na = P.na;
da = P.da;
NM = P.gM/P.muM;
NH = 1;

switch state
    case 'init' %
        SH = 0.97*P.PH_stable*NH; %0.9*NH/na/da*ones(na,1); % cell averages
        EH = 0.01*P.PH_stable*NH; %0.1/na/da*ones(na,1);
        DH = 0.01*P.PH_stable*NH;
        AH = 0.01*P.PH_stable*NH;
        
        % for mosquitoes - assume at equilibrium
        PH = SH+EH+DH+AH;
        [SM,EM,IM] = mosquito_ODE(DH,AH,PH,NM);
        
        Cm = 0*ones(na,1);
        Cac = 0*ones(na,1);
        Ctot = P.c1*Cac+P.c2*Cm;
        
    case 'EE' % start from EE    
        [SH, EH, DH, AH, Cac, Cm, Ctot] = steady_state('EE','numerical');
        PH = SH+EH+DH+AH;
        [SM,EM,IM] = mosquito_ODE(DH,AH,PH,NM);        

    otherwise
        keyboard
end

P.phi = sigmoid_prob(Ctot./P.PH_stable, 'phi'); % prob. of DH -> RH
P.rho = sigmoid_prob(Ctot./P.PH_stable, 'rho'); % prob. of severely infected, EH -> DH
P.psi = sigmoid_prob(Ctot./P.PH_stable, 'psi'); % prob. AH -> DH

end