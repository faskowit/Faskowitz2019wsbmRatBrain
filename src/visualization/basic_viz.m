%% raw data

plotWSBM(templateModel.Data.Raw_Data,ones(154,1)')
axis square
xlabel('')
ylabel('')

%%
subplot(4,2,1:4)

plotWSBM(templateModel)
axis square
xlabel('')
ylabel('')

subplot(4,2,5)
imagesc(reshape(templateModel.Para.predict_e,6,6)')
axis square
xlabel('predicted edge-existence')
colorbar

subplot(4,2,6)
imagesc(reshape(templateModel.Para.predict_w,6,6)')
axis square
xlabel('predicted weight')
colorbar

subplot(4,2,7)
scatter(templateModel.Para.predict_e,templateModel.Para.predict_w)
axis square
cb = colorbar();
cb.Visible = 'off'
xlabel('predicted edge-existence')
ylabel('predicted weight')

subplot(4,2,8)
%function [weiBM,avgWeiBM,binBM,avgBinBM,stdWeiBM,sizesMat] = get_block_mat(CIJ,ca,excludeNaN)
[~,tmp] = get_block_mat(templateAdj,ca_wsbm);
imagesc(tmp)
colorbar
axis square
xlabel('average weight btwn blocks')

%% plot modularity

subplot(1,2,1)

plotWSBM(templateModel)
axis square
xlabel('')
ylabel('')
title('WSBM')

subplot(1,2,2)

ca_mod_match = hungarianMatch(ca_wsbm,ca_mod) ;
plotWSBM(templateModel.Data.Raw_Data,dummyvar(ca_mod_match)')
axis square
xlabel('')
ylabel('')
title('MOD')






