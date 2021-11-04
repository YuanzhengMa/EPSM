% %% 作图1 aline
% D = D3d_3d./max(max(max(D3d_3d)));
% L = size(D,1);
% W = size(D,2);
% H = size(D,3);
% mylegend = {};
% for l = 1:10:200
%     mylegend=['Aline at ',num2str(l),'th slow axis ',num2str(100),'th fast axis'];
%     x = l.*ones(4096,1);
%     y =1:4096;
%     z = squeeze(D3d_sfh(l,100,:));
%     figure(200),plot3(x,y,z,'DisplayName',mylegend,'Color',[rand rand rand]),...
%     legend;
%     hold on
% end
% xlabel('The index of slow axis','Rotation',20,'FontName','Times New Roman');
% ylabel('The points detected along Aline','Rotation',-30,'FontName','Times New Roman');
% zlabel('Amplitude');
%% 作图2 bscan
D = data_3d./max(max(max(data_3d)));
L = size(D,1);
W = size(D,2);
H = size(D,3);
H = 200;
for l = 1:60:L    
    
    [y,z] = meshgrid(1:W,1:H);
    Bscan = squeeze(D(l,1:W,1:H));
    Bzero = zeros(size(Bscan));
    X = l.*ones(size(Bscan,2),l);
    mylegend = ['This is my',num2str(l)]
    x = squeeze(X(:,l));
    figure(201),h = surf(x,y,z,Bscan),h.FaceColor='interp',h.EdgeColor = 'none',h.FaceAlpha=0.9;
    colormap(hot),caxis([0 1]);hold on%'FaceColor','texturemap','FaceAlpha',0.5

end
xlabel('The index of slow axes','Rotation',20,'FontName','Times New Roman');
ylabel('The index of fast axes','Rotation',-30,'FontName','Times New Roman');
zlabel('The index of axis direction','FontName','Times New Roman');
view(gca,[-37.5 30]);
grid(gca,'on');
% 设置其余坐标区属性
set(gca,'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'FontName','Times New Roman');