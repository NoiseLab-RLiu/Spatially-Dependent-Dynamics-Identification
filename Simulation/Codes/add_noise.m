rng(0);
noise = sqrt(2)*randn(size(Waves));
Waves = Waves+noise;
% figure
% imagesc(Waves(:,:,501));
% colorbar
% axis square
% set(gca,'YDir','normal')
% clims = [-30,30];