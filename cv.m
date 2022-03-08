%% Inputs
close all
clear all

nb=load('brain.mat');
I = nb.K;

%I=imread('terminal.gif');
I=im2double(I);

its = 500;
mu = 0.1;
r = 10;A
dt=0.5;
Energy = zeros(its,1);


%-- End
%% Contours
% Bubbles
[M,N] = size(I);
spacing = round(r/4);
basis_box = zeros(2*r+spacing);
basis_circle = zeros(2*r);
basis_quart = zeros(r);
basis_quart(1,ceil(r/2)+1:end) = 1;
basis_quart(2:end,ceil(r/2)+1:end) = 1;
basis_quart(ceil(r/2)+1:end,1) = 1;
col = ceil(r/2);
if r > 2
    for row = 2:ceil(r/2)
        basis_quart(row,col)=1;
        basis_quart((row+1):end,col) = 1;
        col = col-1;

    end
end
basis_circle(1:r,1:r) = basis_quart;
basis_circle(1:r,r+1:end) = rot90(basis_quart,3);
basis_circle(r+1:end,1:r) = rot90(basis_quart);
basis_circle(r+1:end,r+1:end) = rot90(basis_quart,2);
basis_box(1:2*r,1:2*r) = basis_circle;
%basis_box((2*r+spacing):end,(2*r+spacing):end) = 1;
basis_box = basis_box-0.5;
[m,bc] = contour(basis_box);

row_copies = ceil(M/(2*r+spacing));
col_copies = ceil(N/(2*r+spacing));
bubbles = zeros((2*r+spacing)*row_copies,(2*r+spacing)*col_copies);

for rc = 1:row_copies
    for cc = 1:col_copies
        bubbles(((rc-1)*(2*r+spacing)+1):rc*(2*r+spacing),((cc-1)*(2*r+spacing)+1):cc*(2*r+spacing))=basis_box;
    end
end
bubbles = bubbles(1:M,1:N);

% Single 
single_contour = zeros(size(I));
%{
m1 = ceil(1+(floor(0.9*M))*rand(1));
m2 = ceil(m1+(floor(0.9*M)-m1)*rand(1));
n1 = ceil(1+(floor(0.9*N))*rand(1));
n2 = ceil(n1+(floor(0.9*N)-n1)*rand(1));
%single_contour(m1:m2,n1:n2) = rand(1);
single_contour(10,100:170) = 1;
single_contour(50,100:170) = 1;
single_contour(10:50,170) = 1;
single_contour(10:50,100) = 1;
single_contour = single_contour-0.5;
%}
cx = round(M/2);
cy = 50;
R = 50;
[x,y] = meshgrid(1:max(M,N));
circ = (x-cx).^2+(y-cy).^2 <= R^2;
circ = circ-0.5;
single_contour = circ(1:M,1:N);


%{
figure();
subplot(2,2,1);imagesc(basis_quart); 
subplot(2,2,2);imagesc(basis_circle);
subplot(2,2,3);imagesc(basis_box);
subplot(2,2,4);imagesc(bubbles);

figure();
subplot(2,2,1);contour(basis_quart,[0,0]); 
subplot(2,2,2);contour(basis_circle,[0,0]);
subplot(2,2,3);contour(basis_box,[0,0]);
subplot(2,2,4);contour(bubbles,[0,0]);

figure();
subplot(1,2,1);imagesc(single_contour); 
subplot(1,2,2);contour(single_contour,[0,0]);
%}

%-- End
%% Alg
contour_type = 's'; % s (single) / b (bubbles)
if lower(contour_type) == 's'
    phi0 = single_contour;
else
    phi0 = bubbles;
end

% Values at t = 0
F = eps; 
phi = phi0; 

figure();
subplot(1,2,1); imshow(I);% hold on; contour(phi0,'b');
hold on;
contour(phi, [0 0], 'r');
title(['0/' num2str(its) ' Iterations']); 
hold off;
drawnow;

for n=1:its
  % Heaviside and Delta funtions
  % H1
  %{
  H_eps=10^(-5);  
  H = zeros(size(I)); D = zeros(size(I));
  ind= find(phi>H_eps);
  H(ind)=1;
  ind1 = find(phi<H_eps & phi>-H_eps);
  for i=1:length(ind1)
      H(ind1(i))=1/2*(1+phi(ind1(i))/H_eps+1/pi*sin(pi*phi(ind1(i))/H_eps));
      D(ind1(i))=1/(2*H_eps)*(1+cos(pi*phi(ind1(i))/H_eps))+1;
  end
  %D=1/(2*H_eps).*(1+cos(pi/H_eps.*phi));
  %}

  % H2
  %%{
  H_eps=0.1;  
  H = zeros(size(I)); %D = zeros(size(I));
  ind= find(phi>H_eps);
  H(ind)=1;
  ind1 = find(phi<H_eps & phi>-H_eps);
  for i=1:length(ind1)
      H(ind1(i))=1/2*(1+(2/pi)*atan(phi(ind1(i))/H_eps));
  end

 % H = 1/2.*(1+(2/pi).*atan(phi./H_eps));
  D = (H_eps/pi)./(H_eps^2+phi.^2);            
  %%}

  % Curvature
  %phi_pad = padarray(phi,[1,1],1,'pre');
  %phi_pad = padarray(phi_pad,[1,1],1,'post'); 
  phi_pad = padarray(phi,[1,1],1,'both'); % getting the 'ghost' points

  % central difference
  fy = (phi_pad(3:end,2:N+1)-phi_pad(1:M,2:N+1))/2;
  fx = (phi_pad(2:M+1,3:end)-phi_pad(2:M+1,1:N))/2;
  fyy = phi_pad(3:end,2:N+1)+phi_pad(1:M,2:N+1)-2*phi;
  fxx = phi_pad(2:M+1,3:end)+phi_pad(2:M+1,1:N)-2*phi;
  fxy = (1/4).*(phi_pad(3:end,3:end)-phi_pad(1:M,3:end)+phi_pad(3:end,1:N)-phi_pad(1:M,1:N));
  K = ((fxx.*fy.^2-2*fxy.*fx.*fy+fyy.*fx.^2)./((fx.^2+fy.^2+eps).^(3/2))).*(fx.^2+fy.^2).^(1/2);
  K(1,:) = eps;
  K(end,:) = eps;
  K(:,1) = eps;
  K(:,end) = eps;
  K = K./max(max(abs(K)));

  % Evolution
  c1 = sum(sum(I.*H))/(sum(sum(phi>=0))+eps); % equation (6)
  c2 = sum(sum(I.*(1-H)))/(sum(sum(phi<0))+eps); % equation (7) 
  F = (mu*K-(I-c1).^2+(I-c2).^2); % equation (9)
  F = F./max(max(abs(F)));

  phi = phi+dt.*D.*F;
  %phi = phi+dt.*F; %This could work by equation (9) = 0. So we don't really need delta ?? More "noise"
 
  subplot(1,2,1);imshow(I);% hold on; contour(phi0,'b');
  hold on;
  contour(phi, [0 0], 'r');
  title([num2str(n) '/' num2str(its) ' Iterations']); 
  hold off;
  drawnow;

  Energy(n) = norm(D.*F,2);

end
imshow(I); %hold on; contour(phi0,'b');
hold on;
contour(phi, [0 0], 'r');
title([num2str(n) '/' num2str(its) ' Iterations']); 
hold off;

%Segmented image
segment = phi>0; 

subplot(1,2,2); imshow(segment); title('Segmented Image');

figure(); plot(Energy);
