clear all

%http://thevirtualheart.org/FHN2dindex.html

%IMPORTANT. READ BEFORE EXECUTING.

%Constants defined in the line 19 can be generated with functions
%normal_heart(), reentrant_waves() or fibrilation(). Different simulations 
%are obtained from chosen constants.

%normal_heart() simulates a healthy heart.
%reentrant_waves() simulates the reentrant waves.
%fibrilation() simulates a ventricular fibrilations based on reentrant
%waves

%If fibrilation() is used, uncomment the if condition in lines 156-162. 
%For normal_heart() and fibrilation(), those lines must be commented.

[Nt,L,nu,c2,d,t,sigma,alfa,u0,r0,stop] = normal_heart();

N = 50; M = 50; %Use even N and M 
hx = L/(N+1); hy = L/(M+1);
dt=0.1;
theta=1/2; %theta=0 implicit, theta=1 explicit

Beta = @(x,y) dt*sigma(x,y);
x = 0:hx:L; y = 0:hy:L;

%Recovery variable R is initialized with 0. Potential U has an initial
%pulse with value u0 in the top.
R = zeros(N+2, M+2);
U = zeros(N+2, M+2);

for j=1:N+2
    for i=1:M+2
        R(j,i)=0;
        beta(j,i) = Beta(x(j),y(i));
        n = (j-1)*(M+2)+i;
        if j>(N+2)-10
            U(n)=u0;
        end
    end
end

beta = reshape(beta',[(N+2)*(M+2),1]);
U=reshape(U,[N+2,M+2])';

%Neumann Condition
U(1,1) = U(2,2); U(1,M+2) = U(2,M+1);
U(N+2,1) = U(N+1,2); U(N+2,M+2) = U(N+1,M+1);
for i=2:N
   U(i,1)=U(i,2); U(i,M+2)=U(i,M+1); 
end
for j=2:M
   U(1,j)=U(2,j); U(N+2,j)=U(N+1,j); 
end

%% Ax Matrix for first derivative with respect to x 
Ax = 1/hx*gallery('tridiag',zeros(1,(M+2)*(N+2)-1),-1*ones(1,(M+2)*(N+2)),ones(1,(M+2)*(N+2)-1));

%Backward derivative is used on the boundaries
for i=0:N
    Ax(i*(M+2)+M+2,i*(M+2)+(M+2)+1)=0;
    Ax(i*(M+2)+M+2,i*(M+2)+(M+2))=1/hx;
    Ax(i*(M+2)+M+2,i*(M+2)+(M+2)-1)=-1/hx;
end
Ax((N+1)*(M+2)+M+2,(N+1)*(M+2)+(M+2))=1/hx;
Ax((N+1)*(M+2)+M+2,(N+1)*(M+2)+(M+2)-1)=-1/hx;

%First derivative with respect to x of beta
sgm_x = Ax*beta;

%% Axx Matrix for second derivative with respect to x
Axx = (1/hx^2)*gallery('tridiag',ones(1,(M+2)*(N+2)-1),-2*ones(1,(M+2)*(N+2)),ones(1,(M+2)*(N+2)-1));

for i=0:N
   Axx(i*(M+2)+(M+2),i*(M+2)+(M+2)+1)=0;
   Axx(i*(M+2)+(M+2)+1,i*(M+2)+(M+2))=0;
end

Ax2 = speye((M+2)*(N+2));

%% Ay Matrix for first derivative with respect to y

aux = gallery('tridiag',zeros(1,(N+2)-1),-1*ones(1,(N+2)),ones(1,(N+2)-1));
Ay = (1/hy)*kron(aux,speye(M+2)); 

%Backward derivative is used on the boundaries
for j=1:M+2
    Ay((N+1)*(M+2)+j,(N+1)*(M+2)+j) = 1/hy;
    Ay((N+1)*(M+2)+j,(N)*(M+2)+j) = -1/hy;
end
%Se calcula la primera derivada en y de betay
sgm_y = Ay*beta;

%% Ayy Matrix for second derivative with respect to y

aux = gallery('tridiag',ones(1,(N+2)-1),-2*ones(1,(N+2)),ones(1,(N+2)-1));
Ayy = (1/hy^2)*kron(aux,speye(M+2));


%% Final part

B1 = Ax2-(1-theta)*(diag(beta)*Axx+diag(beta)*Ayy+diag(sgm_x)*Ax+diag(sgm_y)*Ay);
B2 = (1-theta)*dt*c2*speye((N+2)*(M+2));
B3 = -((1-theta)*dt*c2)*speye((N+2)*(M+2));
B4 = (1+(1-theta)*dt*c2*d)*speye((N+2)*(M+2));

B = [B1 , B2 ; B3 , B4];

C1 = Ax2+theta*(diag(beta)*Axx+diag(beta)*Ayy+diag(sgm_x)*Ax+diag(sgm_y)*Ay);
C2 = (-theta)*dt*c2*speye((N+2)*(M+2));
C3 = ((theta)*dt*c2)*speye((N+2)*(M+2));
C4 = (1-(theta)*dt*c2*d)*speye((N+2)*(M+2));

C = [C1,C2;C3,C4];

%Plotting initial state
figure('units','normalized','outerposition',[0 0 1 1])
x = 0:hx:L; y = 0:hy:L;
[x,y]=meshgrid(x,y);
subplot(1,2,1)
contourf(x,y,U,'LineColor','None')
zlim([-u0-0.5,u0+0.5])
colorbar
caxis([-u0-0.5,u0+0.5])
title(strcat('Solution for u(x,y,t), N = ',num2str(N),', M = ',num2str(M),' and \Delta t = ',num2str(dt)))
xlabel('x [cm]'); ylabel('y [cm]'); zlabel('u(x,y)')
pause(0.00005)
%Electrocardiogram
ecg = [U((N+2)/2+15,(M+2)/2)-1/2*(U((N+2)/2,M+2-10)+U((N+2)/2,1+10))];
subplot(1,2,2)
xlabel('time [ms]');ylabel('Voltage[mV]')
plot([0],ecg)
ylim([-1,1]); xlim([0,Nt*dt])

for i=1:Nt   
    U = reshape(U',[(N+2)*(M+2),1]);
    R = reshape(R',[(N+2)*(M+2),1]);
 
    %The following if creates a new pulse from top to bottom (do not
    %comment)
    if mod(i,t)==0 && i<stop
        for j = 1:N+2
            for k = 1:M+2
                n = (j-1)*(M+2)+k;
                if j>(N+2)-10
                    U(n)=u0;
                end
            end
        end
    end
    
    %The following lines create pulses from left to right. Uncomment if
    %fibrilation() is used
%     if mod(i,t)==15 && i<stop
%         for j = 0:N+1
%             for k = 1:10
%                 U(j*(M+2)+k)=u0;
%             end
%         end
%     end
    
    %Solving the system for time n+1
    V = B\(C*[U ; R] + [dt.*nu.*U.*(U-alfa).*(1-U); zeros((M+2)*(N+2),1)]); 
    U = V(1:(M+2)*(N+2)); 
    R = V((M+2)*(N+2)+1:end);
    U = reshape(U,[N+2,M+2])';
    
    %Neumann condition
    U(1,1) = U(2,2); U(1,M+2) = U(2,M+1);
    U(N+2,1) = U(N+1,2); U(N+2,M+2) = U(N+1,M+1);
    for l=2:N
       U(l,1)=U(l,2); U(l,M+2)=U(l,M+1); 
    end
    for j=2:M
       U(1,j)=U(2,j); U(N+2,j)=U(N+1,j); 
    end
    
    R = reshape(R,[N+2,M+2])';
    
    %Electrocardiogram
    ecg = [ecg, U((N+2)/2+15,(M+2)/2)-1/3*(U((N+2)/2,M+2-10)+U((N+2)/2,1+10)+U(10,(M+2)/2))];
    
    %Plotting solution U
    subplot(1,2,1)
    xlabel('x [cm]'); ylabel('y [cm]');
    contourf(x,y,U,'LineColor','None');
    title(strcat('Solution for u(x,y,t), N = ',num2str(N),', M = ',num2str(M),' and \Delta t = ',num2str(dt)))
    xlabel('x'); ylabel('y'); zlabel('u(x,y)')
    grid on
    colorbar
    caxis([-u0-0.5,u0+0.5])
    zlim([-u0-0.5,u0+0.5])
    pause(0.005)
    
    %Plotting ECG
    subplot(1,2,2)
    plot(100*dt*(0:i),ecg)
    xlim([0,Nt*dt*100])
    ylim([-1,1])
    xlabel('Time [ms]'); ylabel('Voltage[mV]')
    title('Electrocardiogram')
end