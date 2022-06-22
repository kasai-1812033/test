%% Initialize
clc; clear; close all;
N = 3; %number of agents
dim = 3; %次元
p = 0.2*rand(dim*N,1); %初期位置
p = [1 0 0 0 1 0 0 0 1]';
P = p'; %初期位置の転置（1行1ステップ）
alpha = 1; %正定値

%% calculation
%重み関数
phi  = @(q,p) min(sqrt(sum((q - reshape(p,dim,N)).^2))); %ℓ2ノルム
% phi  = @(q,p) min(sum(abs(q - reshape(p,dim,N)))); %ℓ1ノルム

%重み関数内の引数pを取り除く
%arrayfunはベクトル化を行うMATLAB内の標準関数
if dim==2
    phiv = @(X,Y,p) arrayfun(@(x,y) phi([x,y]',p),X,Y);
else
    phiv = @(X,Y,Z,p) arrayfun(@(x,y,z) phi([x,y,z]',p),X,Y,Z)
end


%%
%評価関数（pを除いたΦを積分）
%quad2dは2次元関数の数値積分を行うMATLAB11内の標準関数
if dim==2
    J = @(p) quad2d(@(X,Y) phiv(X,Y,p),0,1,0,1,'Singular',false); %2次元均一配置分布
else
    J = @(p) integral3(@(X,Y,Z) phiv(X,Y,Z,p),0,0.5,0,0.5,0,0.5); %3次元均一配置分布
end

%入力式のdJdpの計算
%numjacはヤコビアンを数値計算するMATLAB内の標準関数
for i=1:N
    px(i,:)=p(3*i-2); %エージェントのx座標
    py(i,:)=p(3*i-1); %エージェントのy座標
    pz(i,:)=p(3*i); %エージェントのz座標
end
Ps = -[1 1 1]+ 3*[0,0,0;0,1,0;1,0,0;1,1,0;0,0,1;0,1,1;1,0,1;1,1,1]; %ボロノイ分割用の座標ベクトル
Ps = [px,py,pz;Ps];
[v,c] = voronoin(Ps); %3次元ボロノイ分割
for i=1:N
    [k{i},av{i}] = convhull([v(c{i},1),v(c{i},2),v(c{i},3)]); %エージェント周りのボロノイ空間
    V{i} = alphaShape([v(c{i},1),v(c{i},2),v(c{i},3)]); %エージェント周りのボロノイ空間
end
% 質量
mass = @(p,i) integral3(@(X,Y,Z) phiv(X,Y,Z,p), ...
    min(V{i}.Points(:,1)),max(V{i}.Points(:,1)), ...
    min(V{i}.Points(:,2)),max(V{i}.Points(:,2)), ...
    min(V{i}.Points(:,3)),max(V{i}.Points(:,3)));
% 重心
centx = @(p,i) integral3(@(X,Y,Z) X.*phiv(X,Y,Z,p), ...
    min(V{i}.Points(:,1)),max(V{i}.Points(:,1)), ...
    min(V{i}.Points(:,2)),max(V{i}.Points(:,2)), ...
    min(V{i}.Points(:,3)),max(V{i}.Points(:,3)))/mass(p,i);
centy = @(p,i) integral3(@(X,Y,Z) Y.*phiv(X,Y,Z,p), ...
    min(V{i}.Points(:,1)),max(V{i}.Points(:,1)), ...
    min(V{i}.Points(:,2)),max(V{i}.Points(:,2)), ...
    min(V{i}.Points(:,3)),max(V{i}.Points(:,3)))/mass(p,i);
centz = @(p,i) integral3(@(X,Y,Z) Z.*phiv(X,Y,Z,p), ...
    min(V{i}.Points(:,1)),max(V{i}.Points(:,1)), ...
    min(V{i}.Points(:,2)),max(V{i}.Points(:,2)), ...
    min(V{i}.Points(:,3)),max(V{i}.Points(:,3)))/mass(p,i);
% dJdp = @(p) numjac( @(t,p) J(p),0,p,J(p),1e-6*ones(dim*N,1),[],0); %標準関数
% dJdp = @(p) 2*mass(p)*(p-cent(p));

%% plot
if dim==2
    % 2次元用
    for step=1:100
        clf
        fcontour(@(X,Y) phiv(X,Y,p),[0,1,0,1],'r','LevelStep',1e-2); %エージェント周りの等高線
        axis equal
        hold on
        title(sprintf('k=%d',step)); %1ステップごとにタイトルに表示
        plot(P(:,1:2:end-1),P(:,2:2:end),'.-'); %エージェントの位置
        drawnow %この時点で描画
    
        p = p - alpha * dJdp(p)'; %状態更新
        P(end+1,:) = p; %更新した位置の追記
    end
else
    % 3次元用
    for step=1:100
        for i=1:N
            clf
        %     fsurf(@(X,Y,Z) phiv(X,Y,p)); %エージェント周りの等高線
            axis equal
            grid on
            hold on
            view(10,-25);
            title(sprintf('k=%d',step)); %1ステップごとにタイトルに表示
            plot3(P(:,1:3:end-2),P(:,2:3:end-1),P(:,3:3:end),'.-'); %エージェントの位置
            drawnow %この時点で描画

            cent_x(i)= centx(p,i)
            cent_y(i)= centy(p,i)
            cent_z(i)= centz(p,i)
            volume(i)= mass(p,i)
            p(3*i-2:3*i) = p(3*i-2:3*i) - alpha * (2*volume(i)*(p(3*i-2:3*i)-[cent_x(i);cent_y(i);cent_z(i)])); %状態更新でmassとcentを使って計算
            P(end+1,:) = p; %更新した位置の追記
            
            
            plot3(cent_x(i),cent_y(i),cent_z(i),'o');
            if i==1
                trisurf(k{i},v(c{i},1),v(c{i},2),v(c{i},3),'Facecolor',[1,0,0],'Facealpha',alpha*0.3)
            elseif i==2
                trisurf(k{i},v(c{i},1),v(c{i},2),v(c{i},3),'Facecolor',[0,1,0],'Facealpha',alpha*0.3)
            else
                trisurf(k{i},v(c{i},1),v(c{i},2),v(c{i},3),'Facecolor',[0,0,1],'Facealpha',alpha*0.3)
            end
        end
    end
end