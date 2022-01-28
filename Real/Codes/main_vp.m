load('Xsus.mat');

%%

A_raw = zeros(9,8100,5);
C_res = nan(100,100,5);
percent = zeros(5,1);% success rate
me = zeros(5,1);
ma = zeros(5,1);
mp = zeros(5,1);
for cfreq=3:3
    VP_Butter6;
    TwoD_dict_int;
    lasso_seq;
    A_raw(:,:,cfreq-2) = a_raw;
    recover_c_a_2dwave;
    C_res(:,:,cfreq-2) = c_rec;
end
%% Plotting
% for i=1:5
% figure
% imagesc(C_res(:,:,i),'AlphaData',~isnan(C_res(:,:,i)))
% colorbar
% caxis([450 700])
% axis square
% ylabel('y (mm)','interpreter','latex','FontSize',22)
% xlabel('x (mm)','interpreter','latex','FontSize',22)
% title(strcat(num2str((i+2)*10), ' kHz'), 'interpreter','latex','FontSize',22)
% % xticks([10 30 ])
% % xticklabels({'30','40','50','60','70'})
% ylim([5,95])
% xlim([5,95])
% ax=gca
% ax.TickLabelInterpreter = 'latex';
% ax.FontSize=28;
% set(gca, 'YDir','normal')
% 
% % imagesc((((abs(A_raw(:,:,i))))))
% % % cbh = colorbar('XTickLabel',{'0','0.5','1'}, ...
% % %                'XTick',0:0.5:1,'TickLabelInterpreter','latex')
% % colorbar 
% % caxis([0,0.4])
% % ylabel('Index of the entry $i$ in ${\bar{\mathbf{a}}}_n$','interpreter','latex')
% % xlabel('Location index $n$','interpreter','latex')
% % title(strcat('$|{\bar{\mathbf{a}}}_n(i)|$ for', {' '}, num2str(10*(i+2)), ' kHz'),'interpreter','latex')%
% % yticklabels({'1','2','3','4','5','6','7','8','9'})
% % %xticks([1,200,400,600,784])
% % %xticks([1,10,20,30,40,51])
% % %xticks([1,20,40,60,80,97])
% % ax=gca
% % ax.FontSize=20;
% % %caxis([-5 0])
% % colormap(hot)
% % set(gca,'TickLabelInterpreter','latex')
% end

 %% Calculating mean and variation
% meanc = zeros(5,1);
% varc = zeros(5,1);
% madc = zeros(5,1);
% for i=1:5
%     C_tmp = C_res(:,:,i);
%     meanc(i) = mean(C_tmp(:),'omitnan');
%     varc(i) = var(C_tmp(:),'omitnan');
%     C_diff = C_tmp - meanc(i)*ones(100,100);
%     madc(i) = sum(sum(abs(C_diff),'omitnan'),'omitnan')/(sum(sum((~isnan(C_tmp)))));
%     %varc(i) = nanvar(C_tmp(i));
% end

