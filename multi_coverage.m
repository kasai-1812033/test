%% Initialize
clc; clear; close all;
N = 3; %number of agents
dim = 3; %次元
p = 0.2*rand(dim*N,1); %初期位置
p = [1 0 0 0 1 0 0 0 1]';
P = p'; %初期位置の転置（1行1ステップ）
alpha = 1; %正定値

%% calculation
for step=1:100
    if dim == 2
        for i = 1:N
            px(i,:)=p(3*i-2); %エージェントのx座標
    py(i,:)=p(3*i-1); %エージェントのy座標
for i=1:N
    px(i,:)=p(3*i-2); %エージェントのx座標
    py(i,:)=p(3*i-1); %エージェントのy座標
    pz(i,:)=p(3*i); %エージェントのz座標
end
Ps = -[1,1,1]+ 3*[0,0,0;0,1,0;1,0,0;1,1,0;0,0,1;0,1,1;1,0,1;1,1,1]; %ボロノイ分割用の座標ベクトル
Ps = [px,py,pz;Ps];
[v,c] = voronoin(Ps); %3次元ボロノイ分割
for i=1:N
    [k{i},av{i}] = convhull(v(c{i},1),v(c{i},2),v(c{i},3),'Simplify',true); %エージェント周りのボロノイ空間
    TR = triangulation(k{i},[v(c{i},1),v(c{i},2),v(c{i},3)]);
    F = faceNormal(TR);
    Ptri = incenter(TR);
    
    d = 0.05;% gridの刻み幅
    [qx,qy,qz] = meshgrid(-3:d:3,-3:d:3,-3:d:3);
    bx = [reshape(qx,[numel(qx),1]),reshape(qy,[numel(qx),1]),reshape(qz,[numel(qx),1])];
        
    % 質量
%     高速版
    zo = find(max(sum(Ptri.*F,2)-(F*bx')<0,[],1)==0);
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
    cent(i,:) = dmass/mass
    logger.cent{step} = cent;
    
    axis equal
    grid on
    hold on
    view(10,-25);
    title(sprintf('k=%d',step)); %1ステップごとにタイトルに表示
    plot3(P(step,1:3:end-2),P(step,2:3:end-1),P(step,3:3:end),'.-'); %エージェントの位置
    drawnow %この時点で描画
    if i == 1
        trisurf(k{i},v(c{i},1),v(c{i},2),v(c{i},3),'Facecolor','c','Facealpha',alpha*0.01);
    elseif i==2
        trisurf(k{i},v(c{i},1),v(c{i},2),v(c{i},3),'Facecolor','r','Facealpha',alpha*0.01);
    else
        trisurf(k{i},v(c{i},1),v(c{i},2),v(c{i},3),'Facecolor','g','Facealpha',alpha*0.01);
    end
    plot3(cent(i,1),cent(i,2),cent(i,3),'o');

    p(3*i-2:3*i) = p(3*i-2:3*i) - 0.1 * (p(3*i-2:3*i)-[cent(i,1);cent(i,2);cent(i,3)]); %状態更新でmassとcentを使って計算
    P(end+1,:) = p; %更新した位置の追記
end
end