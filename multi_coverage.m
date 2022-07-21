%% Initialize
clc; clear; close all;
N = 3; %number of agents
dim = 3; %次元
% p = 0.2*rand(dim*N,1); %初期位置
p = [1 0 0 0 1 0 0 0 1]';
% p = round(2*rand(1,9))';
P = p'; %初期位置の転置（1行1ステップ）
alpha = 1; %正定値
phi0 = [2 2 2]; %重み位置
ts = 0;
dt = 1/30;
te = 10;
%% calculation
t = ts;
step = 1;
while round(t,5)<=te
    logger.time(step) = t;
    for i = 1:N
        px(i,:)=p(3*i-2); %エージェントのx座標
        py(i,:)=p(3*i-1); %エージェントのy座標
        pz(i,:)=p(3*i); %エージェントのz座標
    end
    Ps = -6*[1,1,1]+ 12*[0,0,0;0,1,0;1,0,0;1,1,0;0,0,1;0,1,1;1,0,1;1,1,1]; %ボロノイ分割用の座標ベクトル
    Ps = [px,py,pz;Ps];
    logger.Ps{step} = Ps;
    [v,c] = voronoin(Ps); %3次元ボロノイ分割
    logger.v{step} = v;
    logger.c{step} = c;
    for i = 1:N
        
        [k{i},av{i}] = convhull(v(c{i},1),v(c{i},2),v(c{i},3),'Simplify',true); %エージェント周りのボロノイ空間
        logger.k{i,step} = k{i};
        logger.av{i,step} = av{i};
        TR = triangulation(k{i},v(c{i},1),v(c{i},2),v(c{i},3)); %三角形分割
        F = faceNormal(TR); %三角形分割した面に対する法線ベクトル
        logger.F{i,step} = F;
        Ptri = incenter(TR); %三角形分割した面の内心
        logger.Ptri{i,step} = Ptri;
        
        d = 0.05;% gridの刻み幅
        [qx,qy,qz] = meshgrid(-2:d:2,-2:d:2,-2:d:2);
        bx = [reshape(qx,[numel(qx),1]),reshape(qy,[numel(qx),1]),reshape(qz,[numel(qx),1])];
        logger.bx{i,step} = bx;
            
        % 質量
    %     高速版
%         zo = find(max(sum(Ptri.*F,2)-(F*bx')<0,[],1)==0);
        input1 = sum(Ptri.*F,2); %三角形分割した面の内心と単位法線ベクトルを掛け合わせたのの合計（ボロノイ分割空間）
        input2 = (F*bx'); %単位法線ベクトルと微小空間を掛け合わせる（ボクセルの位置）
        input3 = max(input1-input2<0,[],1); %ボクセルが重み
        zo = find(input3 == 0); %0になるときだけボクセルが丸ごとボロノイ空間内に存在
        phi_d = normpdf(phi0 - bx(zo,:)); %重み位置と領域内ボクセルとの距離の正規分布関数
        logger.zo{i,step} = zo;
        
%         dmass = sum(bx(zo,:),1);
        dmass = sum(bx(zo,:)'*phi_d,1); %重み付きボクセル
        mass = sum(dmass,2); %全部の重みを合算
%         mass = length(zo);
        
        % 重心
        cent = dmass/mass;
        logger.cent{i,step} = cent;

%         if mod(step,(1/dt)*1) == 1
%             figure(8)
%             title(['t=',num2str(logger.time(step)),'s'],'FontSize',25);
%             hold on
%             
%             view(45,45);
%             daspect([1 1 1]);
%             ax = gca;
%             ax.Box = 'on';
%             ax.GridColor = 'k';
%             ax.GridAlpha = 0.4;
%             xlabel('x [m]','FontSize',12);
%             ylabel('y [m]','FontSize',12);
%             zlabel('z [m]','FontSize',12);
%             xlim([-2,2]);
%             ylim([-2,2]);
%             zlim([-2,2]);
%             grid on
%     
%             plot3(P(step,1:3:end-2),P(step,2:3:end-1),P(step,3:3:end),'k*','MarkerSize',10); %エージェントの位置
%             
%             trisurf(k{i},v(c{i},1),v(c{i},2),v(c{i},3),'Facecolor','r','Facealpha',alpha*0.3);
%     %         plot3(logger.bx{i,time}(logger.zo{i,time},1),logger.bx{i,time}(logger.zo{i,time},2),logger.bx{i,time}(logger.zo{i,time},3),'g.');
%                 
%             
%             plot3(cent(1),cent(2),cent(3),'bo','MarkerSize',10);
%             
%     
%             hold off
%         end
            
        p(3*i-2:3*i) = p(3*i-2:3*i) - 0.1 * (p(3*i-2:3*i)-cent'); %状態更新でmassとcentを使って計算
    end
    P(end+1,:) = p'; %更新した位置の追記
    t = t + dt
    step = step + 1;
end
logger.P = P;
%% 動画
close all;
pt = 15; %FontSize
figure(1)
Fig.No = 1;
time = 1;
v = VideoWriter(['3D voronoiFormation N = ',num2str(N)],'MPEG-4');
v.FrameRate = 30;
open(v);
while time <= numel(logger.time)
    clf(figure(1))
    disp('.');

    hold on
    view(135,35); %normal
%     view(0,90); %top
%     view(0,0); %front
    daspect([1 1 1])
    ax = gca;
    ax.Box = 'on';
    ax.GridColor = 'k';
    ax.GridAlpha = 0.4;
    xlabel('x [m]','FontSize',pt);
    ylabel('y [m]','FontSize',pt);
    zlabel('z [m]','FontSize',pt);
    xlim([-2,2]);
    ylim([-2,2]);
    zlim([-2,2]);
    grid on;

    title(sprintf('t=%fs',logger.time(time)),'FontSize',25); %1ステップごとにタイトルに表示
    plot3(logger.P(time,1:3:end-2),logger.P(time,2:3:end-1),logger.P(time,3:3:end),'k*','MarkerSize',10); %エージェントの位置
    
    
    for i = 1:N
        trisurf(logger.k{i,time},logger.v{time}(logger.c{time}{i},1),logger.v{time}(logger.c{time}{i},2),logger.v{time}(logger.c{time}{i},3),'EdgeColor','none','Facecolor','r','Facealpha',alpha*0.3);
%         plot3(logger.bx{i,time}(logger.zo{i,time},1),logger.bx{i,time}(logger.zo{i,time},2),logger.bx{i,time}(logger.zo{i,time},3),'g.');
            
        
        plot3(logger.cent{i,time}(1),logger.cent{i,time}(2),logger.cent{i,time}(3),'bo','MarkerSize',10);
    end

    hold off

    if mod(time,(1/dt)*2) == 1
        Fig.No = Fig. No + 1;
        figure(Fig.No)
%         title(['t=',num2str(logger.time(time)),'s'],'FontSize',25);
        hold on
        
        view(135,35);
        daspect([1 1 1]);
        ax = gca;
        ax.XAxis.Visible = 'off';
        ax.YAxis.Visible = 'off';
        ax.ZAxis.Visible = 'off';
        ax.Box = 'on';
        ax.GridColor = 'k';
        ax.GridAlpha = 0.4;
%         xlabel('x [m]','FontSize',pt);
%         ylabel('y [m]','FontSize',pt);
%         zlabel('z [m]','FontSize',pt);
        xlim([-2,2]);
        ylim([-2,2]);
        zlim([-2,2]);
        grid on

        plot3(logger.P(time,1:3:end-2),logger.P(time,2:3:end-1),logger.P(time,3:3:end),'k*','MarkerSize',10); %エージェントの位置
        for i = 1:N
            trisurf(logger.k{i,time},logger.v{time}(logger.c{time}{i},1),logger.v{time}(logger.c{time}{i},2),logger.v{time}(logger.c{time}{i},3),'EdgeColor','none','Facecolor','r','Facealpha',alpha*0.3);
    %         plot3(logger.bx{i,time}(logger.zo{i,time},1),logger.bx{i,time}(logger.zo{i,time},2),logger.bx{i,time}(logger.zo{i,time},3),'g.');
                
            
            plot3(logger.cent{i,time}(1),logger.cent{i,time}(2),logger.cent{i,time}(3),'bo','MarkerSize',10);
        end

        hold off
    end
    time = time + 1;

    frame = getframe(figure(1));
    writeVideo(v,frame);
end
disp('simulation end');
close(v);
%% 自己位置と重心の推移
% figure(9) %agent1
% log_cent = cell2mat(logger.cent');
% hold on
% axis equal
% xlabel('t [s]','FontSize',pt);
% ylabel('Position [m]','FontSize',pt)
% xlim([0,logger.time(end)]);
% ylim([-2,2]);
% plot(logger.time,P(1:end-1,1));
% plot(logger.time,P(1:end-1,2));
% plot(logger.time,P(1:end-1,3));
% plot(logger.time,log_cent(:,1),'--');
% plot(logger.time,log_cent(:,2),'--');
% plot(logger.time,log_cent(:,3),'--');
% hold off