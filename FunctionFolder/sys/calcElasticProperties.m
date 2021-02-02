function elasticProperties=calcElasticProperties(elasticMatrix, numIons, atomType, Volume)

%%%%%%%   This funciton is used to calculate the elastic tensor related properites of the structure
%%%%%%%   
%%%%%%%   Added by Haiyang Niu and Guangrui Qian
		

%% Units 
	h = 6.62606957e-34; %% m2*kg/s    Planck`s constant
	k = 1.3806488e-23;  %% J/K        Boltzmann`s constant
	NA=	6.0221413e+23;  %% 1/mol      Avogadro`s number 
 


%%% When we have a empty elestic constant tensor, we have such result:
if isempty(elasticMatrix)
	elasticProperties=[NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
	return;
end

%%% Calcuate the structure density	
density = calcDensity(numIons, atomType, Volume);	
	
% C refers to elastic stiffness tensor for easy typo
C = elasticMatrix;

% Equations to calculate 
 % * Kv  bulk modulus Kv, 
 % * Gv	 shear modulus Gv, 
 % * Ev  Young`s modulus , 
 % * vv  Poisson`s ratio v using the Voigt method
Kv = ((C(1,1)+C(2,2)+C(3,3)) + 2*(C(1,2)+C(2,3)+C(3,1)))/9;
Gv = ((C(1,1)+C(2,2)+C(3,3)) - (C(1,2)+C(2,3)+C(3,1)) + 3*(C(4,4)+C(5,5)+C(6,6)))/15;
Ev = 1/((1/(3*Gv))+(1/(9*Kv)));
vv  = 0.5*(1-((3*Gv)/(3*Kv+Gv)));

S = inv(C);      %S refers to elastic compliance tensor, S is the inverse matrix of C
% Equations to calculate bulk modulus Kv, shear modulus Gv, Young`s modulus Ev, Poisson`s ratio v using the Reuss method
Kr = 1/((S(1,1)+S(2,2)+S(3,3))+2*(S(1,2)+S(2,3)+S(3,1)));
Gr = 15/(4*(S(1,1)+S(2,2)+S(3,3))-4*(S(1,2)+S(2,3)+S(3,1))+3*(S(4,4)+S(5,5)+S(6,6)));
Er = 1/((1/(3*Gr))+(1/(9*Kr)));
vr = 0.5*(1-((3*Gr)/(3*Kr+Gr)));

% Equations to calculate bulk modulus Kv, shear modulus Gv, Young`s modulus Ev, Poisson`s ratio v and Vicker`s hardness
K_h = (Kv+Kr)/2;   % Bulk modulus
G_h = (Gv+Gr)/2;   % Shear modulus
E_h = (Ev+Er)/2;   % Young`s modulus
v_h = (vv+vr)/2;   % Poisson`s ratio
GK  = G_h/K_h;     % Pugh`s modulus ratio, shear modulus over bulk modulus

%% Mechanical stability condition
% eighvalues of the elastic stiffness tensor matrix much be positive
E = eig(C) ; % eighvalues of the elastic stiffness tensor matrix much be
is_stable = all(E>0);  


%% constraint for the Bulk and Shear modulus
%if G_h < 1e-3 
%   G_h = 0;
%end
%if K_h < 1e-3
%   K_h = 0;	   
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
if G_h*K_h < 1e-3
  Hv = NaN;
else
  Hv = 2*((G_h*G_h*G_h)/(K_h*K_h))^0.585-3;   % Equation to calculate Vicker`s hardness by Chen-Niu hardness model
end
	  
% Calculating longitudinal and transverse elastic wave velocity from Navier`s equation
V_p = ((K_h+(4*G_h/3))/density)^0.5;         % Velocity S-wave
V_s = (G_h/density)^0.5;                     % Velocity P-wave
V_m = (1/3*(2/(V_s^3)+1/(V_p^3)))^(-1/3);    % Mean sound velocity


%Calculating the Debye temperature       
atomicMass=zeros(1, length(atomType));
for i = 1:length(atomType)
	atomicMass(i) = elementMass(atomType(i));
end	  
    n = sum(numIons);                               % n is the number of atoms in the structure
Sum_M = sum(numIons.*atomicMass);                   % Sum_M is the summanry of atomic weight
 De_T = h/k*((3*n)/(4*pi)*(NA*density/Sum_M))^(1/3)*V_m*100000;    % Debye temperature

%  Calculating the Fracture toughness of 3D solid material
Kg =0 ;

%                   1    2    3    4   5   6   7    8     9    10   11    12
elasticProperties=[K_h, G_h, E_h, v_h, GK, Hv, Kg, De_T, V_m, V_s, V_p, is_stable];
