% LOCALIZAÇÃO DO PONTO DE FALHA POR MODELO DE IMPEDÂNCIA
% Equacionamento para valores de Vx (tensão no ponto de falta na blindagem)
% Simulação de falta a terra na Fase A
% Cálculo para valores instantâneos 

% Parâmetros da blindagem metálica dos cabos isolados:

% Rb = 3.146e-3 * 1000; % Resistência total da blindagem
% Lb = 9.44e-7 * 1000; % Indutância total da blindagem
% C = 0.36 * 1e-6; % Capacitância total do cabo

% Dados de entrada manuais para o circuito elétrico

% valor_usuario = input('Digite o valor da Resistência da Blindagem Metálica em [ohm/m]: ');
% valor_usuario = input('Digite o valor da Indutância da Blindagem Metálica em [H/m]: ');
% valor_usuario = input('Digite o valor da Capacitância do cabo [F/km]: ');

% Vetores para armazenar os resultados de x:

n = length(tout);
x1 = zeros(1, n-2); 
x2 = zeros(1, n-2);

% Calculando as derivadas numéricas:

% Fundamental

dt = diff(tout); % Intervalos de tempo
dVLA_VbA_dt = diff(VLA - VbA) ./ dt; % Derivada de (VLA - VbA) em relação ao tempo
dILA_dt = diff(ILA) ./ dt; % Derivada de ILA em relação ao tempo
dVsA_dt = diff(VsA) ./ dt; % Derivada de VsA em relação ao tempo
dIsA_dt = diff(IsA) ./ dt; % Derivada de IsA em relação ao tempo
d2VsA_dt2 = diff(dVsA_dt) ./ dt(1:end-1); % Segunda derivada de VsA
d2VLA_VbA_dt2 = diff(dVLA_VbA_dt) ./ dt(1:end-1); % Segunda derivada de (VLA - VbA)
dIbA_dt = diff(IbA) ./dt; % Derivada de IbA em relação ao tempo

% Sequencia Zero

dVL0_Vb0_dt = diff(VL0 - Vb0) ./ dt; % Derivada de (VL0 - Vb0) em relação ao tempo
dIL0_dt = diff(IL0) ./ dt; % Derivada de IL0 em relação ao tempo
dVs0_dt = diff(Vs0) ./ dt; % Derivada de Vs0 em relação ao tempo
dIs0_dt = diff(Is0) ./ dt; % Derivada de Is0 em relação ao tempo
d2Vs0_dt2 = diff(dVs0_dt) ./ dt(1:end-1); % Segunda derivada de Vs0
d2VL0_Vb0_dt2 = diff(dVL0_Vb0_dt) ./ dt(1:end-1); % Segunda derivada de (VL0 - Vb0)
dIb0_dt = diff(Ib0) ./dt; % Derivada de Ib0 em relação ao tempo

% Loop para calcular x para cada ponto de dados:

for i = 1:n-2
    
    % Coeficientes da equação quadrática para x (fundamental):
    
    A = - (Rb*C/2)*(dVsA_dt(i)) - (Lb*C/2)*(d2VsA_dt2(i)) - (Rb*C/2)*(dVLA_VbA_dt(i)) - (Lb*C/2)*(d2VLA_VbA_dt2(i));
    B = Rb*IbA(i) + Lb*dIbA_dt(i) + Rb*C*(dVLA_VbA_dt(i)) + Lb*C*(d2VLA_VbA_dt2(i));
    C_term = - VbA(i) - (Rb*C/2)*(dVLA_VbA_dt(i)) - (Lb*C/2)*(d2VLA_VbA_dt2(i));

    % Raízes de x com valores instantâneos:

    if A ~= 0
        x1(i) = real((-B + sqrt(B.^2 - 4 * A * C_term)) / (2 * A));
        x2(i) = real((-B - sqrt(B.^2 - 4 * A * C_term)) / (2 * A));
    else
        x1(i) = -C_term/B;
        x2(i) = NaN;
    end

    % Coeficientes da equação quadrática para x na Sequencia Zero:
    
    A0 = - (Rb*C/2)*(dVs0_dt(i)) - (Lb*C/2)*(d2Vs0_dt2(i)) - (Rb*C/2)*(dVL0_Vb0_dt(i)) - (Lb*C/2)*(d2VL0_Vb0_dt2(i));
    B0 = Rb*Ib0(i) + Lb*dIb0_dt(i) + Rb*C*(dVL0_Vb0_dt(i)) + Lb*C*(d2VL0_Vb0_dt2(i));
    C0_term = - Vb0(i) - (Rb*C/2)*(dVL0_Vb0_dt(i)) - (Lb*C/2)*(d2VL0_Vb0_dt2(i));

    % Raízes de x com valores sequencia zero:

    if A0 ~= 0
        x10(i) = real((-B0 + sqrt(B0.^2 - 4 * A0 * C0_term)) / (2 * A0));
        x20(i) = real((-B0 - sqrt(B0.^2 - 4 * A0 * C0_term)) / (2 * A0));
    else
        x10(i) = -C0_term/B0;
        x20(i) = NaN;
    end
end

% Gráfico dos valores das raízes da equação (x):

figure;
subplot(2,1,1)
plot(tout(1:end-2), abs(x1), '-o', ... % Adiciona marcadores 'o' (círculos)
     tout(1:end-2), abs(x2), '-s', 'LineWidth', 1); % Adiciona marcadores 's' (quadrados)
xlabel('Time (s)');
ylabel('(pu)');
% axis([0.004199 0.004203 0 1])
title('Fundamental Values Location');
grid on;

subplot(2,1,2)
plot(tout(1:end-2), abs(x10), '-o', ... % Adiciona marcadores 'o' (círculos)
     tout(1:end-2), abs(x20), '-s', 'LineWidth', 1); % Adiciona marcadores 's' (quadrados)
xlabel('Time (s)');
ylabel('(pu)');
% axis([0.004199 0.004203 0 1]);
title('Zero Sequence Location');
grid on;

