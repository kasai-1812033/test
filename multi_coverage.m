%% Initialize
clc; clear; close all;
N = 8; %number of agents
dim = 3; %次元
p = 0.1*rand(dim*N,1); %初期位置
P = p'; %初期位置の転置（1行1ステップ）
alpha = 1; %正定値
%% voronoi
% clc; clear; close all;
% x=[0.2,0.6,0.9,0.7,0.8,0.5,0.4,0.1,0.3,0.5];
% y=[0.7,0.9,0.8,0.6,0.3,0.4,0.2,0.1,0.3,0.5];
% x = P(:,1:2:end);
% y = P(:,2:2:end);
% X = [x',y'];
% [VX,VY] = voronoi(x,y);
% h = plot(VX,VY,'-b',X(:,1),X(:,2),'ro');
% 
% % Assign labels to the points.
% nump = size(X,1);
% plabels1 = arrayfun(@(n) {sprintf('x%d', n)}, (1:nump)');
% hold on
% Hp1 = text(X(:,1)-0.05, X(:,2), plabels1, 'Color','r','FontWeight', ...
%       'bold','HorizontalAlignment','center', ...
%       'BackgroundColor', 'none');
% 
% % Assign labels to the voronoi areas.
% plabels2 = arrayfun(@(n) {sprintf('C%d', n)}, (1:nump)');
% Hp2 = text(X(:,1)+0.03, X(:,2)-0.03, plabels2, 'Color','b','FontWeight', ...
%       'bold', 'HorizontalAlignment','center', ...
%       'BackgroundColor', 'none');
% hold off
% 
% axis equal
% xlim([0,1]);
% ylim([0,1]);
%% calculation
%重み関数
phi  = @(q,p) min(sqrt(sum((q - reshape(p,dim,N)).^2))); %ℓ2ノルム
% phi  = @(q,p) min(sum(abs(q - reshape(p,dim,N)))); %ℓ1ノルム

%重み関数内の引数pを取り除く
%arrayfunはベクトル化を行うMATLAB内の標準関数
if dim==2
    phiv = @(X,Y,p) arrayfun(@(x,y) phi([x,y]',p),X,Y);
else
    phiv = @(X,Y,Z,p) arrayfun(@(x,y,z) phi([x,y,z]',p),X,Y,Z);
end

%評価関数（pを除いたΦを積分）
%quad2dは2次元関数の数値積分を行うMATLAB11内の標準関数
if dim==2
    J = @(p) quad2d(@(X,Y) phiv(X,Y,p),0,1,0,1,'Singular',false); %2次元均一配置分布
else
    J = @(p) integral3(@(X,Y,Z) phiv(X,Y,Z,p),0,0.5,0,0.5,0,0.5); %3次元均一配置分布
end

%入力式のdJdpの計算
%numjacはヤコビアンを数値計算するMATLAB内の標準関数
dJdp = @(p) numjac( @(t,p) J(p),0,p,J(p),1e-6*ones(dim*N,1),[],0);

%% plot
if dim==2
    % 2次元用
    for k=1:100
        clf
        fcontour(@(X,Y) phiv(X,Y,p),[0,1,0,1],'r','LevelStep',1e-2); %エージェント周りの等高線
        axis equal
        hold on
        title(sprintf('k=%d',k)); %1ステップごとにタイトルに表示
        plot(P(:,1:2:end-1),P(:,2:2:end),'.-'); %エージェントの位置
        drawnow %この時点で描画
    
        p = p - alpha * dJdp(p)'; %状態更新
        P(end+1,:) = p; %更新した位置の追記
    end
else
    % 3次元用
    for k=1:100
        clf
    %     fsurf(@(X,Y,Z) phiv(X,Y,p)); %エージェント周りの等高線
        axis equal
        grid on
        hold on
        view(10,-25);
        title(sprintf('k=%d',k)); %1ステップごとにタイトルに表示
        plot3(P(:,1:3:end-2),P(:,2:3:end-1),P(:,3:3:end),'.-'); %エージェントの位置
        drawnow %この時点で描画
    
        p = p - alpha * dJdp(p)'; %状態更新
        P(end+1,:) = p; %更新した位置の追記
    end
end