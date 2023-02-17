% This function gives the corresponding force vector given the displacement vector using the Bouc-Wen class model proposed in the following reference
% Reference: Oh, S., Kim, T., & Song, J. (2023). Bouc–Wen class models considering hysteresis mechanism of RC columns in nonlinear dynamic analysis.

% The parameters below are defined under the scale factor A of the original Bouc-Wen model is an unit (A=1).
%{

k0         :   Normalized initial stiffness
Fy         :   Normalized yield force
alpha      :   Post-yield stiffness ratio
beta       :   Basic hysteresis shape control parameter
gamma      :   Basic hysteresis shape control parameter (gamma = 1-beta is used. See the reference for more information.)
n          :   Softening parameter. Controls the transition from linear to non-linear range
deltaNu    :   Strength degradation rate. deltaNu = 0 means no strength degradation.
deltaEta   :   Stiffness degradation rate. deltaEta = 0 means no stiffness degradation.
pin_zeta0  :   Pinching parameter about the total slip
pin_p      :   Pinching parameter about the pinching slope
pin_q      :   Pinching parameter about the location of the pinching initiation
pin_psi    :   Pinching parameter about the pinching magnitude
pin_delpsi :   Pinching parameter about the pinching rate
pin_lambda :   Pinching parameter about the pinching severity
c_eps      :   Deterioration energy amplification factor
c_pch      :   Crack closure coefficient

%}

function force = BoucWen(parameters, displacement)

% Model parameter settings
k0 = (2*pi/parameters(1))^2/9.8; % kN/mm or g/m
Fy = parameters(2); % kN or g
alpha = parameters(3);
beta = parameters(4);
gamma = 1-beta;
n = parameters(5);
deltaNu = parameters(6);
deltaEta = parameters(7);
pin_zeta0 = parameters(8);
pin_p = parameters(9);
pin_q = parameters(10);
pin_psi = parameters(11);
pin_delpsi = parameters(12);
pin_lambda = parameters(13);
c_eps = parameters(14);
c_h = parameters(15);

dispYield = Fy/k0;

% Basic parameters for the algorithm
maxIter = 50000;        % maximum number of iterations for the Newton-Raphson algorithm
tolerance = 1e-10;      % prescribed tolerance for the convergence check
startingPoint = 1e-5;   % for numerical stability
lr = 0.9;               % learning rate

% Start
force = zeros(length(displacement),1);
displacement = [0; displacement];

DispT = 0.0;
z_prev = 0.0;
e_post_old = 0.0;
e_nega_old = 0.0;
dispPeak = 0;
deltaDisp = displacement(1)-0;

sympref('HeavisideAtOrigin',0);

for ii = 2:length(displacement)

    DispTdT = displacement(ii);
    dir_pri = sign(deltaDisp);
    deltaDisp = DispTdT - DispT;
    dir_post = sign(deltaDisp);
    flag_peak = (dir_pri~=dir_post);
    dispPeak = flag_peak*DispTdT + (1-flag_peak)*dispPeak;

    % Difference between the current displacement and the previous maximum/minimum displacement
    post_del = (displacement(ii)-max(displacement(1:(ii-1))));
    nega_del = (-displacement(ii)+min(displacement(1:(ii-1))));

    % Hysteretic energy amplification factor
    A_post = 1+c_eps*heaviside(post_del);
    A_nega = 1+c_eps*heaviside(nega_del);


    % Perform Newton-Rhapson
    count = 0;
    count_total = 0;
    z_new = 0.0;
    z_old = startingPoint;
    z_eval = startingPoint;

    while abs(z_old - z_new) > tolerance && count < maxIter

        % Step 1 - Calculate the evaluation function Gamma
        e_post_new = e_post_old + A_post*(1 - alpha)*deltaDisp*k0/Fy*z_eval;
        e_nega_new = e_nega_old + A_nega*(1 - alpha)*deltaDisp*k0/Fy*z_eval;
        e_new = heaviside(displacement(ii))*e_post_new + (1-heaviside(displacement(ii)))*e_nega_new;

        nu_new = 1 + deltaNu*e_new;
        eta_new = 1 + deltaEta*e_new;

        if nu_new < 0
            lr = lr*0.1;
            count_total = count_total+count;
            count = 0;
            z_eval = startingPoint;

            e_post_new = e_post_old + A_post*(1 - alpha)*deltaDisp*k0/Fy*z_eval;
            e_nega_new = e_nega_old + A_nega*(1 - alpha)*deltaDisp*k0/Fy*z_eval;
            e_new = heaviside(displacement(ii))*e_post_new + (1-heaviside(displacement(ii)))*e_nega_new;

            nu_new = 1 + deltaNu*e_new;
            eta_new = 1 + deltaEta*e_new;
        end

        a_1 = beta*sign(deltaDisp*z_eval)+ gamma;
        a_2 = (1-abs(z_eval)^n*a_1*nu_new)/eta_new;

        Zu = (1/((beta+gamma)*nu_new))^(1/n);
        zeta1 = pin_zeta0*(1-exp(-pin_p*e_new))*(1-exp(-c_h*(abs(dispPeak)/dispYield)));
        zeta2 = (pin_psi+e_new*pin_delpsi)*(pin_lambda+zeta1);

        h = 1 - zeta1*exp(-(z_eval*sign(deltaDisp)-pin_q*Zu)^2/zeta2^2);

        Gamma = z_eval - z_prev - h*a_2*deltaDisp*k0/Fy;

        % Step 2 - Evaluate the deriviative of Gamma with respect to z
        e_post_new_ = A_post*(1 - alpha)*k0/Fy*deltaDisp;
        e_nega_new_ = A_nega*(1 - alpha)*k0/Fy*deltaDisp;
        e_new_ = heaviside(displacement(ii))*e_post_new_ + (1-heaviside(displacement(ii)))*e_nega_new_;

        nu_new_ = deltaNu*e_new_;
        eta_new_ = deltaEta*e_new_;
        Zu_ = -nu_new_*(beta+gamma)/n*((beta+gamma)*nu_new)^(-(n+1)/n);

        a_2_ = (-eta_new_*(1-abs(z_eval)^(n)*a_1*nu_new) ...
            -eta_new*(n*abs(z_eval)^(n-1)*a_1*nu_new*sign(z_eval) ...
            +abs(z_eval)^(n)*a_1*nu_new_))/eta_new^2;

        zeta1_ = pin_zeta0*pin_p*exp(-pin_p*e_new)*e_new_*(1-exp(-c_h*(abs(dispPeak)/dispYield)));
        zeta2_ = pin_psi*zeta1_ + pin_lambda*pin_delpsi*e_new_ + pin_delpsi*e_new_*zeta1 + pin_delpsi*e_new*zeta1_;

        a3 = -exp(-(z_eval*sign(deltaDisp)-pin_q*Zu)^2/zeta2^2);
        a4 = 2*zeta1*(z_eval*sign(deltaDisp)-pin_q*Zu)*(sign(deltaDisp)-pin_q*Zu_)/zeta2^2;
        a5 = 2*zeta1*(z_eval*sign(deltaDisp)-pin_q*Zu)^2/zeta2^3;
        h1_ = a3*(zeta1_ - a4 + zeta2_*a5);

        Gamma_ = 1 - (h1_*a_2+h*a_2_)*deltaDisp*k0/Fy;

        % Step 3 - Obtain the trial value of z
        z_new = z_eval - lr*Gamma/Gamma_;

        % Step 4 - Update the trial value
        z_old = z_eval;
        z_eval = z_new;

        count = count + 1;

        % Warning if there is no convergence
        if count == maxIter
            disp(['WARNING: Convergence failed after max iteration. Learning rate will decrease by 10%. Current learning rate is ',num2str(lr)]);
            lr = lr*0.1;
            count = 0;
            z_new = 0.0; z_eval = startingPoint;
        end
    end

    % Compute resisting force.
    fs = alpha*k0*DispTdT + (1 - alpha)*Fy*z_eval;
    force(ii-1) =  fs;

    DispT = DispTdT;
    z_prev = z_eval;
    e_post_old = e_post_new;
    e_nega_old = e_nega_new;
end
end
