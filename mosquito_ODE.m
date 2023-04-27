function [SM,EM,IM] = mosquito_ODE(DH,AH,PH,NM)
% evolve mosquito quantities, keep at the steady state 
global P

[~,bM] = biting_rate(PH,NM);
NH = trapz(PH)*P.da;
lamM = FOI_M(bM,DH,AH,NH);
SM = P.gM/(lamM+P.muM);
IM = (P.gM/P.muM)*(P.sigma/(P.sigma+P.muM))*(lamM/(lamM+P.muM));
EM = (P.muM/P.sigma)*IM;

end