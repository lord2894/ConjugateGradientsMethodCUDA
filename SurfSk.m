N = 1000;
M = importdata('C:\\Users\\Kir\\Documents\\outputGrid.txt')
ZCGM=zeros(N);
Y = M(1:N,2);
%Читаем данные
for i=1:N
    ZCGM(i,:) = M(N*(i-1)+1:N*i,3);
    X(i) = M(N*(i-1)+1,1);
end
X = X.';
%Строим график численного решения
figure
surfCGM = surf(X,Y,ZCGM);
set(surfCGM,'LineStyle','none')
set(surfCGM,'FaceColor',[0 0 1],'FaceAlpha',0.5);
xlabel('X') 
ylabel('Y')
zlabel('P')
grid on
% Вычисляем точное решение
phi = @(x,y)exp(1-x*x*y*y);
for idx = 1:numel(X)
    x = X(idx);
    for idy = 1:numel(Y)
        y=Y(idy);
        ZPHI(idx,idy) = phi(x,y);
    end
end
% Строим график точного решения
figure
surfPHI = surf(X,Y,ZPHI);
set(surfPHI,'LineStyle','none')
set(surfPHI,'FaceColor',[1 0 0],'FaceAlpha',0.5);
grid on
xlabel('X') 
ylabel('Y')
zlabel('U')
%Строим графики точного и численного решения
figure
surfCGM2 = surf(X,Y,ZCGM);
set(surfCGM2,'LineStyle','none')
set(surfCGM2,'FaceColor',[0 0 1],'FaceAlpha',0.5);
grid on
hold on
surfPHI2 = surf(X,Y,ZPHI);
set(surfPHI2,'LineStyle','none')
set(surfPHI2,'FaceColor',[1 0 0],'FaceAlpha',0.5);
xlabel('X') 
ylabel('Y')
zlabel('U|P')
%Вычисляем невязку
ZDiscrepancy = abs(ZPHI - ZCGM);
ZDiscrepancyV = reshape(ZDiscrepancy',N*N,1)';
% Строим график невязки
figure
%C = gradient(ZDiscrepancy);
surfDiscrepancy = surf(X,Y,ZDiscrepancy);
set(surfDiscrepancy,'LineStyle',':')
n = norm(ZDiscrepancyV)
xlabel('X') 
ylabel('Y')
zlabel('Psi')
