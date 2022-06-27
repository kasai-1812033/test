[x1,y1,z1] = sphere(24);
x1 = x1(:);
y1 = y1(:);
z1 = z1(:);
x2 = x1+5;
P = [x1 y1 z1; x2 y1 z1];
P = unique(P,'rows');
plot3(P(:,1),P(:,2),P(:,3),'.')
axis equal
grid on

shp = alphaShape(P(:,1),P(:,2),P(:,3),1);
plot(shp)
axis equal
xp=rand(10,1);
yp = rand(10,1);
zp = rand(10,1);
[tf,ID] = inShape(shp,xp,yp,zp);
%%
close all
[x,y,z] = meshgrid(0:1:2,-2:1:0,-1:1:1);
x = x(:);
y = y(:);
z = z(:);
x = x1+0.1;
y = y1+0.2;
z = z1+0.3;

[k1,av1] = convhull(x,y,z);

% trisurf(k1,x,y,z,'FaceColor','cyan')
% axis equal

[k2,av2] = convhull(x,y,z,'Simplify',true);

trisurf(k2,x,y,z,'FaceColor','cyan')
axis equal

TR = triangulation(k2,[x,y,z]);
V = vertexNormal(TR);
F = faceNormal(TR);
P = incenter(TR);
%%
close all
trisurf(TR,'FaceColor',[0.8 0.8 1.0],'FaceAlpha',0.3);
axis equal
hold on
% quiver3(x,y,z, ...
%      V(:,1),V(:,2),V(:,3),0.5,'Color','b');


%  quiver3(P(:,1),P(:,2),P(:,3), ...
%       F(:,1),F(:,2),F(:,3),0.5,'color','r');

%p=6*rand(10,3)-3;
dd = 0.1;
dx = dd;
dy = dd;
dz = dd;
[px,py,pz] = meshgrid([-3:dx:3],[-3:dy:3],[-3:dz:3]);
p = [reshape(px,[numel(px),1]),reshape(py,[numel(px),1]),reshape(pz,[numel(px),1])];
% plot3(p(:,1),p(:,2),p(:,3),'ro')
%%
N = 0;
M = 0;
j = 0;
tic
for i = 1:length(p)
    if isempty(find(sum((P-p(i,:)).*F,2) < 0,1))
        plot3(p(i,1),p(i,2),p(i,3),'b*')
        dm = dx*dy*dz; % *rho(p(i,1),p(i,2),p(i,3));
        N = N + dm*[p(i,1),p(i,2),p(i,3)];
        M = M + dm;
        j = j+1;
    end
end
j
COG = N/M
toc
%% 高速化バージョン
tic
zo=find(max(sum(P.*F,2)-(F*p')<0,[],1)==0);
plot3(p(zo,1),p(zo,2),p(zo,3),'ro')
N = sum(p(zo,:),1);
M = length(zo);
COG = N/M
toc