function [SM,EM,IM] = mosquito_ODE(DH,AH,NH,NM)
% evolve mosquito quantities, keep at the steady state 
global P

[~,bM] = biting_rate(NH,NM);
lamM = FOI_M(bM,DH,AH);
SM = P.gM/(lamM+P.muM);
IM = (P.gM/P.muM)*(P.sigma/(P.sigma+P.muM))*(lamM/(lamM+P.muM));
EM = (P.muM/P.sigma)*IM;

end