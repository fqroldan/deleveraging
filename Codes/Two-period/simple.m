%% Code for the speed of deleveraging project
% Francisco Roldán, NYU. 7/2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%                                    %%%%%%%%%%
%%%%%%%%%%   OBSOLETE. REFER TO MAIN_CODE.M   %%%%%%%%%%
%%%%%%%%%%                                    %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameters
global beta beta_borr chi phi aL B0g B0h B1bar G1 G2 Nbar R0S0 w1 w2
beta = (.97)^5; beta_borr = .7^5;
chi = .65;
gamma = 1;
phi = 1;
aL = .4;
w1 = 1; w2 = 1;

B0g = 0.1;
B0h = 0.1;

B1bar = .1;

G1 = 0;
G2 = 0;

Nbar = 1;

R0S0 = B0g/(1-chi) + (chi/(1-chi))*B0h;

% Known functions
U = @(c1,c2,n1,b) log(c1) + b*log(c2) - aL* n1^(1+phi)/(1+phi) - b*aL*Nbar^(1+phi)/(1+phi);
RequiredTransfers = @(T1,R1) -R1*(B0g + T1 + G1) - G2;                  %T2

% Guess from pre-running
%       [ c1b   ,c2b, c1s   , c2s  , s1  , n1    , b1g ]
guess = [ 1.0951, .9, 1.3808,1.1857,.1857, 1.1951,B1bar];
guessb = guess;

% Prepare variables
T1 = linspace(-B0g,0,101);
R1 = linspace(1,1/beta*1.1,75);


C1s = NaN(length(T1),length(R1)); C1b = C1s; C2s = C1b; C2b = C2s;
N1 = C1s; S1 = C1s;
Us = C1s; Ub = C1b; i = 0;

% Fill the blanks
for jt = 1:length(T1)
    t1 = T1(jt);
    s1 = RequiredSavings(t1);
    for jr = 1:length(R1)
        if jr==1
            guess = guessb;
        end
        r1 = R1(jr);
        t2 = RequiredTransfers(t1,r1);
        % compute assuming borrower at constraint
        [ c1b,c2b,c1s,c2s,s1,n1,b1g ] = EulerSaver( t1,r1,t2,guess );
        % if the constraint doesn't bind
        if c1b*beta_borr*r1> c2b
            [ c1b,c2b,c1s,c2s,s1,n1,b1g ] = bothEuler( t1,r1,t2,guess );
        end
        % Save results
        C1b(jt,jr) = c1b;
        C1s(jt,jr) = c1s;
        C2b(jt,jr) = c2b;
        C2s(jt,jr) = c2s;
        S1(jt,jr) = s1;
        N1(jt,jr) = n1;
        % welfare
        Us(jt,jr) = U(C1s(jt,jr),C2s(jt,jr),N1(jt,jr),beta);
        Ub(jt,jr) = U(C1b(jt,jr),C2b(jt,jr),N1(jt,jr),beta_borr);
        guess = [ c1b,c2b,c1s,c2s,s1,n1,b1g ];
        if jr==1
            guessb = guess;
        end
    end
end

%% Some plots
figure(1)
subplot(3,2,1), plot(R1,Us(1,:),'b-',R1,Us(end,:),'r-'), legend('high repayment','low repayment'), title('Welfare of Savers'), xlabel('R_1')
subplot(3,2,2), plot(R1,Ub(1,:),'b-',R1,Ub(end,:),'r-'), legend('high repayment','low repayment'), title('Welfare of Borrowers'), xlabel('R_1')
subplot(3,2,3), plot(R1,S1(1,:),'b--',...
                     R1,C1s(1,:),'b-',...
                     R1,C1b(1,:),'r-',...
                     R1,N1(1,:),'k--'), legend('S1','C_1^s','C_1^b','gdp'), xlabel('R_1'), title('High repayment')
subplot(3,2,4), plot(R1,S1(end,:),'b--',...
                     R1,C1s(end,:),'b-',...
                     R1,C1b(end,:),'r-',...
                     R1,N1(end,:),'k--'), legend('S1','C_1^s','C_1^b','gdp'), xlabel('R_1'), title('Low repayment')
subplot(3,2,5), plot(R1,C2s(1,:),'b-',...
                     R1,C2b(1,:),'r-',...
                     R1,Nbar*ones(size(R1))), legend('C_2^s','C_2^b','gdp'), xlabel('R_1'), title('High repayment')
subplot(3,2,6), plot(R1,C2s(end,:),'b-',...
                     R1,C2b(end,:),'r-',...
                     R1,Nbar*ones(size(R1))), legend('C_2^s','C_2^b','gdp'), xlabel('R_1'), title('Low repayment')
                 
figure(2)
subplot(3,2,1), plot(T1,Us(:,1),'b-',T1,Us(:,end),'r-'), legend('low R_1','high R_1'), title('Welfare of Savers'), xlabel('T_1')
subplot(3,2,2), plot(T1,Ub(:,1),'b-',T1,Ub(:,end),'r-'), legend('low R_1','high R_1'), title('Welfare of Borrowers'), xlabel('T_1')
subplot(3,2,3), plot(T1,S1(:,1),'b--',...
                    T1,C1s(:,1),'b-',...
                    T1,C1b(:,1),'r-',...
                    T1,N1(:,1),'k--'), legend('S1','C_1^s','C_1^b','gdp'), xlabel('T_1'), title('Low R_1')
subplot(3,2,4), plot(T1,S1(:,end),'b--',...
                    T1,C1s(:,end),'b-',...
                    T1,C1b(:,end),'r-',...
                    T1,N1(:,end),'k--'), legend('S1','C_1^s','C_1^b','gdp'), xlabel('T_1'), title('High R_1')
subplot(3,2,5), plot(T1,C2s(:,1),'b-',...
                    T1,C2b(:,1),'r-',...
                    T1,Nbar*ones(size(T1))), legend('C_2^s','C_2^b','gdp'), xlabel('T_1'), title('Low R_1')
subplot(3,2,6), plot(T1,C2s(:,end),'b-',...
                    T1,C2b(:,end),'r-',...
                    T1,Nbar*ones(size(T1))), legend('C_2^s','C_2^b','gdp'), xlabel('T_1'), title('High R_1')