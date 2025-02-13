
a=load('tau-tau1-tau2-DeltaC1-DeltaC2-DeltaC3.txt');

figure(1)
plot(a(:,1),a(:,5),"k-o",'LineWidth', 1, 'MarkerFaceColor', 'k', 'MarkerSize', 10)
hold on
xlabel('S_{1}','FontSize',27)
ylabel('\DeltaCC','FontSize',27);
title('N1-N2')
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
% plot([0.057 0.057],[65.1799868254079,109.1029868254079],'r--','LineWidth',1)
% axis([-0.000092498812546,0.065792001187454,65.1799868254079,109.1029868254079])
ax = gca();
% ax.XRuler.Exponent = -2;
ax.YRuler.Exponent = -2;
print(h, '-r600', '-dpdf', ['DeltaCC_N1-N2','.pdf']);


