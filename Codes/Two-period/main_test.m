%% Code for the optimal speed of deleveraging project
% Francisco Roldán, NYU. 7/2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%                                 %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%    FOR TESTING PURPOSES ONLY    %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%                                 %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all; clear; clc
disp('>> On the Optimal Speed of Deleveraging')
%% Parameters
global beta beta_borr chi phi aL B0g B0h B1bar G1 G2 Nbar w1 Tb Ts
Delev_parameters

%% Find first best
% Guess from pre-running
guess = [-0.0141, 0.0306, -0.0568, 1.0000, 1.0000, 0.4311,1/beta,B1bar/beta];
[ z ] = fsolve(@(x) first_best(x), guess,optimoptions('fsolve','Display','off'));
aL = z(5);
% T = z(1); Tb = z(2); Ts = z(3); C = z(4); S = z(6);
% R = z(7); Bh = z(8);
% 
% clear T Tb Ts Cb Cs S R Bh
disp('   Steady-state first best')
disp('         T        Tb        Ts         C        aL         S         R        Bh')
disp(z)

gdp = 1;
C = z(4);
S = z(6);
Ts = z(3);
Tb = z(2);
clear z
%% Equilibrium for each T1 and R1
% Known functions
RequiredSavings = @(T1) (chi/(1-chi))* B1bar + (B0g + T1 + G1)/(1-chi); %S1
RequiredTransfers = @(T1,R1) -R1*(B0g + T1 + G1) - G2;                  %T2
U = @(c1,c2,n1,b) log(c1) + b*log(c2) - aL* n1^(1+phi)/(1+phi) - b*aL*Nbar^(1+phi)/(1+phi);

% Guess from pre-running (for BothEuler.m)
%       [ c1b   ,c2b, c1s   , c2s  , s1  , n1    , b1h ]
guess = [ 1.0951, .9, 1.3808,1.1857,.1857, 1.1951,B1bar];
guessb = guess;

T1 = linspace(-B0g,0,301);
% T1 = linspace(-0.16,-.14,301);
R1 = linspace(1,1/beta*1.1,25);

C1s = NaN(length(T1),length(R1)); C1b = C1s; C2s = C1b; C2b = C2s;
N1 = C1s; S1 = C1s; C1bE = C1b;
Us = C1s; Ub = C1b;
%% Loop
eulerstop = 0; mktstop = 0;
for jt = 1:length(T1)
    t1 = T1(jt);
    for jr = 1:length(R1)
        if jr==1
            guess = guessb;
        end
        r1 = R1(jr);
        t2 = RequiredTransfers(t1,r1);
        % if the constraint doesn't bind
        [ c1b,c2b,c1s,c2s,s1,n1,b1h ] = bothEuler( t1,r1,t2,guess );
        C1bE(jt,jr) = c1b;
        % if the constraint binds
        if b1h > B1bar*r1
            s1 = RequiredSavings(t1);
            n1 = G1 + (Nbar-G2)/(beta*r1) + (chi/(1-chi))*...
               ((1+1/beta)* (B1bar+B0g+t1+G1) - (B0h+B0g)) + Ts * ((1/(beta*r1))-1);
            c1s = n1 + R0S0 - s1 + t1 + Ts;
            c2s = beta*r1*c1s;
            c2b = Nbar/chi - (1-chi)/chi * c2s - G2/chi;
            temp = Nbar - B1bar*r1 + t2 + Tb; if abs(temp-c2b) > 10^-10, disp('warning'); end
            b1h = B1bar*r1;
            c1b = w1*n1 + b1h/r1 - B0h + t1 + Tb;
            temp = n1/chi - (1-chi)*c1s/chi - G1/chi; if abs(c1b - temp) > 10^-10, disp('ojo'); end
            C1bE(jt,jr) = NaN;
        end
        if c1b*beta_borr*r1 >= c2b
            if b1h > B1bar*r1
                stop
            end
            if ~eulerstop && jr == 1
                disp('Euler')
                disp(['T1 = ',num2str(t1),' R1 = ',num2str(r1)])
                eulerstop = 1;
            end
        end
        C1b(jt,jr) = c1b;
        C1s(jt,jr) = c1s;
        C2b(jt,jr) = c2b;
        C2s(jt,jr) = c2s;
        S1(jt,jr) = s1;
        N1(jt,jr) = n1;
        if 1 < 4*chi*(1-chi)*c1s*temp*aL
            if ~mktstop
                disp('no clearing wages')
                disp(['T1 = ',num2str(t1),' R1 = ',num2str(r1)])
                mktstop = 1;
            end
        end
        % welfare
        Us(jt,jr) = U(C1s(jt,jr),C2s(jt,jr),N1(jt,jr),beta);
        Ub(jt,jr) = U(C1b(jt,jr),C2b(jt,jr),N1(jt,jr),beta);
        guess = [ c1b,c2b,c1s,c2s,s1,n1,b1h ];
        if jr==1
            guessb = guess;
        end
    end
end
%% Prepare output to show

C1bshow = log(C1b) - log(C);
C1bshowE= log(C1bE)- log(C);
C2bshow = log(C2b) - log(C);
C1sshow = log(C1s) - log(C);
C2sshow = log(C2s) - log(C);
S1show = log(S1) - log(S);
Nbarshow = log(Nbar)-log(gdp);
N1show  = log(N1);
%% Some plots
figure(1)
set(1,'position', [0, 150, 1000, 750])
subplot(3,2,1), plot(R1,Us(1,:),'b-',R1,Us(end,:),'r-'), legend('high repayment','low repayment'), title('Welfare of Savers'), xlabel('R_1')
subplot(3,2,2), plot(R1,Ub(1,:),'b-',R1,Ub(end,:),'r-'), legend('high repayment','low repayment'), title('Welfare of Borrowers'), xlabel('R_1')
subplot(3,2,3), hold on
    plot(R1,S1(1,:),'b--')
    plot(R1,C1sshow(1,:),'b-')
    plot(R1,C1bshow(1,:),'r-')
    plot(R1,C1bshowE(1,:),'.','Color',[.5 0 0],'Linewidth',2)
    plot(R1,N1show(1,:),'k--'), legend('S1','C_1^s','C_1^b','C_1^b*','gdp'), xlabel('R_1'), title('High repayment')
subplot(3,2,4), hold on
    plot(R1,S1(end,:),'b--')
    plot(R1,C1sshow(end,:),'b-')
    plot(R1,C1bshow(end,:),'r-')
    plot(R1,C1bshowE(end,:),'.','Color',[.5 0 0],'Linewidth',2)
    plot(R1,N1show(end,:),'k--'), legend('S1','C_1^s','C_1^b','C_1^b*','gdp'), xlabel('R_1'), title('Low repayment')
subplot(3,2,5), plot(R1,C2sshow(1,:),'b-',...
                     R1,C2bshow(1,:),'r-',...
                     R1,Nbarshow*ones(size(R1))), legend('C_2^s','C_2^b','gdp'), xlabel('R_1'), title('High repayment')
subplot(3,2,6), plot(R1,C2sshow(end,:),'b-',...
                     R1,C2bshow(end,:),'r-',...
                     R1,Nbarshow*ones(size(R1))), legend('C_2^s','C_2^b','gdp'), xlabel('R_1'), title('Low repayment')

figure(2)
set(2,'position', [950, 150, 1000, 750])
subplot(3,2,1), plot(T1,Us(:,1),'b-',T1,Us(:,end),'r-'), legend('low R_1','high R_1'), title('Welfare of Savers'), xlabel('T_1')
subplot(3,2,2), plot(T1,Ub(:,1),'b-',T1,Ub(:,end),'r-'), legend('low R_1','high R_1'), title('Welfare of Borrowers'), xlabel('T_1')
subplot(3,2,3), hold on
    plot(T1,S1(:,1),'b--')
    plot(T1,C1sshow(:,1),'b-')
    plot(T1,C1bshow(:,1),'r-')
    plot(T1,C1bshowE(:,1),'.','Color',[.5 0 0],'Linewidth',2)
    plot(T1,N1show(:,1),'k--'), legend('S1','C_1^s','C_1^b','C_1^b*','gdp'), xlabel('T_1'), title('Low R_1')
subplot(3,2,4), hold on
    plot(T1,S1(:,end),'b--')
    plot(T1,C1sshow(:,end),'b-')
    plot(T1,C1bshow(:,end),'r-')
    plot(T1,C1bshowE(:,end),'.','Color',[.5 0 0],'Linewidth',2)
    plot(T1,N1show(:,end),'k--'), legend('S1','C_1^s','C_1^b','C_1^b*','gdp'), xlabel('T_1'), title('High R_1')
subplot(3,2,5), plot(T1,C2sshow(:,1),'b-',...
                    T1,C2bshow(:,1),'r-',...
                    T1,Nbarshow*ones(size(T1))), legend('C_2^s','C_2^b','gdp'), xlabel('T_1'), title('Low R_1')
subplot(3,2,6), plot(T1,C2sshow(:,end),'b-',...
                    T1,C2bshow(:,end),'r-',...
                    T1,Nbarshow*ones(size(T1))), legend('C_2^s','C_2^b','gdp'), xlabel('T_1'), title('High R_1')