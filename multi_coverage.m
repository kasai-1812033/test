%% Initialize
clc; clear; close all;
N = 3; %number of agents
dim = 3; %次元
p = 0.2*rand(dim*N,1); %初期位置
p = [1 0 0 0 1 0 0 0 1]';
P = p'; %初期位置の転置（1行1ステップ）
alpha = 1; %正定値
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
    Ps = -[1,1,1]+ 3*[0,0,0;0,1,0;1,0,0;1,1,0;0,0,1;0,1,1;1,0,1;1,1,1]; %ボロノイ分割用の座標ベクトル
    Ps = [px,py,pz;Ps];
    logger.Ps{step} = Ps;
    [v,c] = voronoin(Ps); %3次元ボロノイ分割
    logger.v{step} = v;
    logger.c{step} = c;
    for i = 1:N
        
        [k{i},av{i}] = convhull(v(c{i},1),v(c{i},2),v(c{i},3),'Simplify',true); %エージェント周りのボロノイ空間
        logger.k{i,step} = k{i};
        logger.av{i,step} = av{i};
        TR = triangulation(k{i},[v(c{i},1),v(c{i},2),v(c{i},3)]);
        F = faceNormal(TR);
        logger.F{i,step} = F;
        Ptri = incenter(TR);
        logger.Ptri{i,step} = Ptri;
        
        d = 0.05;% gridの刻み幅
        [qx,qy,qz] = meshgrid(-3:d:3,-3:d:3,-3:d:3);
        bx = [reshape(qx,[numel(qx),1]),reshape(qy,[numel(qx),1]),reshape(qz,[numel(qx),1])];
        logger.bx{i,step} = bx;
            
        % 質量
    %     高速版
        zo = find(max(sum(Ptri.*F,2)-(F*bx')<0,[],1)==0);
%         input1 = sum(Ptri.*F,2);
%         input2 = (F*bx');
%         input3 = max(input1 - input2 < 0 ,[],1);
%         zo = find(input3 == 0);
        
        dmass = sum(bx(zo,:),1);
        mass = length(zo);
    %     普通版
    %     dmass = 0;
    %     mass = 0;
    %     for j = 1:length(bx)
    %         if isempty(find(sum((Ptri-bx(j,:)).*F,2)<0,1))
    %             dv = d^3;
    %             dmass = dmass + dv*[bx(j,1),bx(j,2),bx(j,3)];
    %             mass = mass + dv;
    %         end
    %     end
    %     log.mass(i) = mass;
        
        % 重心
        cent = dmass/mass;
        logger.cent{i,step} = cent;
            
        p(3*i-2:3*i) = p(3*i-2:3*i) - 0.1 * (p(3*i-2:3*i)-cent'); %状態更新でmassとcentを使って計算
        P(end+1,:) = p; %更新した位置の追記
    end
    t = t + dt
    step = step + 1;
end
logger.P = P;
%% 動画
figure(1)
time = 1;
v = VideoWriter(['3D voronoiFormationSample N = ',num2str(N)],'MPEG-4');
v.FrameRate = 30;
open(v);
while time < numel(logger.time)
    clf(figure(1))
    disp('.');

    hold on
    view(-10,25);
    daspect([1 1 1])
    ax = gca;
    ax.Box = 'on';
    ax.GridColor = 'k';
    ax.GridAlpha = 0.4;
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    xlim([-3,3]);
    ylim([-3,3]);
    zlim([-3,3]);
    grid on;

    title(sprintf('k=%f',logger.time(time))); %1ステップごとにタイトルに表示
    plot3(logger.P(time,1:3:end-2),logger.P(time,2:3:end-1),logger.P(time,3:3:end),'.-'); %エージェントの位置
    
    for i = 1:N
        if i == 1
            trisurf(logger.k{i,time},logger.v{time}(logger.c{time}{i},1),logger.v{time}(logger.c{time}{i},2),logger.v{time}(logger.c{time}{i},3),'Facecolor','c','Facealpha',alpha*0.3);
        elseif i==2
            trisurf(logger.k{i,time},logger.v{time}(logger.c{time}{i},1),logger.v{time}(logger.c{time}{i},2),logger.v{time}(logger.c{time}{i},3),'Facecolor','r','Facealpha',alpha*0.3);
        else
            trisurf(logger.k{i,time},logger.v{time}(logger.c{time}{i},1),logger.v{time}(logger.c{time}{i},2),logger.v{time}(logger.c{time}{i},3),'Facecolor','g','Facealpha',alpha*0.3);
        end
        plot3(logger.cent{i,time}(1),logger.cent{i,time}(2),logger.cent{i,time}(3),'o');
    end

    hold off
    time = time + 1;

    frame = getframe(figure(1));
    writeVideo(v,frame);
end
disp('simulation end');
close(v);