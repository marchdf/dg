%Philip Johnson

clc
clear all

DATA = csvread('Convergence_8_Slim.csv');
eps = 10^(-16);
eJ_best = 1.0;

cJ = 6; %column for Global error

chc = 5; %column for characteristic mesh width

for j = 1:6 
    %standard+upwind
    hc_sU_p1(j) = DATA(j + 0*6 + 0*30, chc);
    eJ_sU_p1(j) = DATA(j + 0*6 + 0*30, cJ);
    hc_sU_p2(j) = DATA(j + 1*6 + 0*30, chc);
    eJ_sU_p2(j) = DATA(j + 1*6 + 0*30, cJ);
    hc_sU_p3(j) = DATA(j + 2*6 + 0*30, chc);
    eJ_sU_p3(j) = DATA(j + 2*6 + 0*30, cJ);
    hc_sU_p4(j) = DATA(j + 3*6 + 0*30, chc);
    eJ_sU_p4(j) = DATA(j + 3*6 + 0*30, cJ);
    hc_sU_p5(j) = DATA(j + 4*6 + 0*30, chc);
    eJ_sU_p5(j) = DATA(j + 4*6 + 0*30, cJ);
    
    %icb+upwind
    hc_iU_p1(j) = DATA(j + 0*6 + 1*30, chc);
    eJ_iU_p1(j) = DATA(j + 0*6 + 1*30, cJ);
    hc_iU_p2(j) = DATA(j + 1*6 + 1*30, chc);
    eJ_iU_p2(j) = DATA(j + 1*6 + 1*30, cJ);
    hc_iU_p3(j) = DATA(j + 2*6 + 1*30, chc);
    eJ_iU_p3(j) = DATA(j + 2*6 + 1*30, cJ);
    hc_iU_p4(j) = DATA(j + 3*6 + 1*30, chc);
    eJ_iU_p4(j) = DATA(j + 3*6 + 1*30, cJ);
    hc_iU_p5(j) = DATA(j + 4*6 + 1*30, chc);
    eJ_iU_p5(j) = DATA(j + 4*6 + 1*30, cJ);
    
    %standard+central
    hc_sC_p1(j) = DATA(j + 0*6 + 2*30, chc);
    eJ_sC_p1(j) = DATA(j + 0*6 + 2*30, cJ);
    hc_sC_p2(j) = DATA(j + 1*6 + 2*30, chc);
    eJ_sC_p2(j) = DATA(j + 1*6 + 2*30, cJ);
    hc_sC_p3(j) = DATA(j + 2*6 + 2*30, chc);
    eJ_sC_p3(j) = DATA(j + 2*6 + 2*30, cJ);
    hc_sC_p4(j) = DATA(j + 3*6 + 2*30, chc);
    eJ_sC_p4(j) = DATA(j + 3*6 + 2*30, cJ);
    hc_sC_p5(j) = DATA(j + 4*6 + 2*30, chc);
    eJ_sC_p5(j) = DATA(j + 4*6 + 2*30, cJ);
    
    %icb+central
    hc_iC_p1(j) = DATA(j + 0*6 + 3*30, chc);
    eJ_iC_p1(j) = DATA(j + 0*6 + 3*30, cJ);
    hc_iC_p2(j) = DATA(j + 1*6 + 3*30, chc);
    eJ_iC_p2(j) = DATA(j + 1*6 + 3*30, cJ);
    hc_iC_p3(j) = DATA(j + 2*6 + 3*30, chc);
    eJ_iC_p3(j) = DATA(j + 2*6 + 3*30, cJ);
    hc_iC_p4(j) = DATA(j + 3*6 + 3*30, chc);
    eJ_iC_p4(j) = DATA(j + 3*6 + 3*30, cJ);
    hc_iC_p5(j) = DATA(j + 4*6 + 3*30, chc);
    eJ_iC_p5(j) = DATA(j + 4*6 + 3*30, cJ);
    
end

% eJ_min = 0.5*eJ_best;
% eJ_min = 10^(-9);

markfacecol={'r';'m';'k';'b'};
marktype={'s';'p';'d';'>'};

figure(1)
% loglog(hc_sU_p1, eJ_trail_p1,'-.','color',[0.5,0.5,0.5],'LineWidth',3); hold on
% loglog(hchar_vs4_p2, eJ_trail_p2,'-.','color',[0.5,0.5,0.5],'LineWidth',3); hold on
% loglog(hchar_vs4_p3, eJ_trail_p3,'-.','color',[0.5,0.5,0.5],'LineWidth',3); hold on
F1 = loglog(hc_sU_p1, eJ_sU_p1,'sr',hc_iU_p1, eJ_iU_p1, 'pm',hc_sC_p1, eJ_sC_p1,'kd',hc_iC_p1,eJ_iC_p1,'>b','MarkerSize',14,'LineWidth',4); 
set(F1,{'Marker'},marktype,{'MarkerFaceColor'},markfacecol);
hold on
xlabel('$\tilde{h}$','Interpreter','LaTex')
ylabel('$E_{G}$','Interpreter','Latex')
title('Global Error, $p=1$','Interpreter','Latex')
legend(F1,'Standard+Upwind','ICBN+Upwind','Standard+Central','ICBN+Central','location','se')
% axis([0.002 0.1 eJ_min 1])
set(gca,'fontsize',26)
axis square
hold off

figure(2)
% loglog(hc_sU_p2, eJ_trail_p2,'-.','color',[0.5,0.5,0.5],'LineWidth',3); hold on
% loglog(hchar_vs4_p2, eJ_trail_p2,'-.','color',[0.5,0.5,0.5],'LineWidth',3); hold on
% loglog(hchar_vs4_p3, eJ_trail_p3,'-.','color',[0.5,0.5,0.5],'LineWidth',3); hold on
F2 = loglog(hc_sU_p2, eJ_sU_p2,'sr',hc_iU_p2, eJ_iU_p2, 'pm',hc_sC_p2, eJ_sC_p2,'kd',hc_iC_p2,eJ_iC_p2,'>b','MarkerSize',14,'LineWidth',4); 
set(F2,{'Marker'},marktype,{'MarkerFaceColor'},markfacecol);
hold on
xlabel('$\tilde{h}$','Interpreter','LaTex')
ylabel('$E_{G}$','Interpreter','Latex')
title('Global Error, $p=2$','Interpreter','Latex')
legend(F2,'Standard+Upwind','ICBN+Upwind','Standard+Central','ICBN+Central','location','se')
% axis([0.002 0.1 eJ_min 1])
set(gca,'fontsize',26)
axis square
hold off

figure(3)
% loglog(hc_sU_p3, eJ_trail_p3,'-.','color',[0.5,0.5,0.5],'LineWidth',3); hold on
% loglog(hchar_vs4_p3, eJ_trail_p3,'-.','color',[0.5,0.5,0.5],'LineWidth',3); hold on
% loglog(hchar_vs4_p3, eJ_trail_p3,'-.','color',[0.5,0.5,0.5],'LineWidth',3); hold on
F3 = loglog(hc_sU_p3, eJ_sU_p3,'sr',hc_iU_p3, eJ_iU_p3, 'pm',hc_sC_p3, eJ_sC_p3,'kd',hc_iC_p3,eJ_iC_p3,'>b','MarkerSize',14,'LineWidth',4); 
set(F3,{'Marker'},marktype,{'MarkerFaceColor'},markfacecol);
hold on
xlabel('$\tilde{h}$','Interpreter','LaTex')
ylabel('$E_{G}$','Interpreter','Latex')
title('Global Error, $p=3$','Interpreter','Latex')
legend(F3,'Standard+Upwind','ICBN+Upwind','Standard+Central','ICBN+Central','location','se')
% axis([0.002 0.1 eJ_min 1])
set(gca,'fontsize',26)
axis square
hold off

figure(4)
% loglog(hc_sU_p4, eJ_trail_p4,'-.','color',[0.5,0.5,0.5],'LineWidth',3); hold on
% loglog(hchar_vs4_p4, eJ_trail_p4,'-.','color',[0.5,0.5,0.5],'LineWidth',3); hold on
% loglog(hchar_vs4_p4, eJ_trail_p4,'-.','color',[0.5,0.5,0.5],'LineWidth',3); hold on
F4 = loglog(hc_sU_p4, eJ_sU_p4,'sr',hc_iU_p4, eJ_iU_p4, 'pm',hc_sC_p4, eJ_sC_p4,'kd',hc_iC_p4,eJ_iC_p4,'>b','MarkerSize',14,'LineWidth',4); 
set(F4,{'Marker'},marktype,{'MarkerFaceColor'},markfacecol);
hold on
xlabel('$\tilde{h}$','Interpreter','LaTex')
ylabel('$E_{G}$','Interpreter','Latex')
title('Global Error, $p=4$','Interpreter','Latex')
legend(F4,'Standard+Upwind','ICBN+Upwind','Standard+Central','ICBN+Central','location','se')
% axis([0.002 0.1 eJ_min 1])
set(gca,'fontsize',26)
axis square
hold off

figure(5)
% loglog(hc_sU_p5, eJ_trail_p5,'-.','color',[0.5,0.5,0.5],'LineWidth',3); hold on
% loglog(hchar_vs4_p5, eJ_trail_p5,'-.','color',[0.5,0.5,0.5],'LineWidth',3); hold on
% loglog(hchar_vs4_p5, eJ_trail_p5,'-.','color',[0.5,0.5,0.5],'LineWidth',3); hold on
F5 = loglog(hc_sU_p5, eJ_sU_p5,'sr',hc_iU_p5, eJ_iU_p5, 'pm',hc_sC_p5, eJ_sC_p5,'kd',hc_iC_p5,eJ_iC_p5,'>b','MarkerSize',14,'LineWidth',4); 
set(F5,{'Marker'},marktype,{'MarkerFaceColor'},markfacecol);
hold on
xlabel('$\tilde{h}$','Interpreter','LaTex')
ylabel('$E_{G}$','Interpreter','Latex')
title('Global Error, $p=5$','Interpreter','Latex')
legend(F5,'Standard+Upwind','ICBN+Upwind','Standard+Central','ICBN+Central','location','se')
% axis([0.002 0.1 eJ_min 1])
set(gca,'fontsize',26)
axis square
hold off


