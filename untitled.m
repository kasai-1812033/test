[x,y,z] = meshgrid(-2:1:2,-2:1:2,-2:1:2);
x = x(:);
y = y(:);
z = z(:);

[k1,av1] = convhull(x,y,z);

trisurf(k1,x,y,z,'FaceColor','cyan')
axis equal